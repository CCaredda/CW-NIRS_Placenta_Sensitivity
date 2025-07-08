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
import matplotlib.gridspec as gridspec

main_path = "/home/caredda/Videos/PROSPEKT/"


class info_subject:
    def __init__(self, ID = np.array([]),
                    Gestation = np.array([]),
                    thickness_skin = np.array([]),
                    thickness_adipose = np.array([]),
                    thickness_muscle=np.array([]),
                    thickness_placenta=np.array([]),
                    distance_to_placenta=np.array([]),
                    total_thickness=np.array([]),
                    Fitzpatrick_scale=np.array([]),
                    ):
        self.ID = ID
        self.Gestation = Gestation
        self.thickness_skin = thickness_skin
        self.thickness_adipose = thickness_adipose
        self.thickness_muscle = thickness_muscle
        self.thickness_placenta = thickness_placenta
        self.distance_to_placenta = distance_to_placenta
        self.total_thickness = total_thickness
        self.Fitzpatrick_scale = Fitzpatrick_scale


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



def read_thickness_data(data_patient,col=''):

    #Load file
    file = main_path + "Placental Segmentation measurements.xlsx"
    data = pd.read_excel(file, skiprows=1)

    # Accessing column names
    column_names = data.columns

    # Get subject id
    subject_id = data['Study Number'].values#[4:])#.astype(int)
    #Gestation
    Gestation = np.asarray(data['Gestation'+col].values,dtype='<U16')
    #thickness skin
    thickness_skin = data['A (SKIN)'+col].values
    #thickness adipose
    thickness_adipose = data['B (ADIPOSE)'+col].values
    #thickness muscle
    thickness_muscle = data['C (MUSCLE)'+col].values
    #thickness placenta
    thickness_placenta = data['D (PLACENTA)'+col].values
    #Distance to placenta
    Distance_to_placenta = data['E (SKIN TO TOP OF PLACENTA)'+col].values
    #total thickness
    total_thickness = data['F (SKIN TO BOTTOM OF PLACENTA)'+col].values

    idx_to_save = np.where(Gestation != 'nan')[0]

    #Convert subject id PMS00ID in int( ID)
    subject_id = subject_id[idx_to_save]
    for i in range (subject_id.shape[0]):
        subject_id[i] = subject_id[i][4:]
    subject_id = subject_id.astype(int)

    #Only get the gestation week
    Gestation = Gestation[idx_to_save]
    for i in range (Gestation.shape[0]):
        Gestation[i] = Gestation[i][0:2]


    thickness_skin = thickness_skin[idx_to_save]
    thickness_adipose = thickness_adipose[idx_to_save]
    thickness_muscle = thickness_muscle[idx_to_save]
    thickness_placenta = thickness_placenta[idx_to_save]
    Distance_to_placenta = Distance_to_placenta[idx_to_save]
    total_thickness = total_thickness[idx_to_save]



    data_patient.ID = np.append(data_patient.ID,subject_id)
    data_patient.Gestation = np.append(data_patient.Gestation,Gestation)
    data_patient.thickness_skin = np.append(data_patient.thickness_skin,thickness_skin)
    data_patient.thickness_adipose = np.append(data_patient.thickness_adipose, thickness_adipose)
    data_patient.thickness_muscle = np.append(data_patient.thickness_muscle, thickness_muscle)
    data_patient.thickness_placenta = np.append(data_patient.thickness_placenta, thickness_placenta)
    data_patient.distance_to_placenta = np.append(data_patient.distance_to_placenta, Distance_to_placenta)
    data_patient.total_thickness = np.append(data_patient.total_thickness, total_thickness)


    # Load skin tones
    file = main_path + "skin tones.xlsx"
    data = pd.read_excel(file)

    # Accessing column names
    column_names = data.columns
    ID = data['Study ID'].values
    Fitzpatrick_scale = data['Fitzpatrick score'].values

    data_patient.Fitzpatrick_scale = np.zeros(data_patient.ID.shape)

    for i in range(ID.shape[0]):
        id = np.where(ID[i] == data_patient.ID)[0]
        if id.size == 0:
            continue
        data_patient.Fitzpatrick_scale[id] = Fitzpatrick_scale[i]



    #Remove nan values


    return data_patient

def convert_f_mel_Fitzpatrick_scale(f_mel):
    scale = ""
    if f_mel == 0.0255:
        scale = "I-II"
    if f_mel == 0.155:
        scale = "III-V"
    if f_mel == 0.305:
        scale = "VI"
    return scale


def convert_Fitzpatrick_scale_f_mel(scale):
    f_mel = np.array([])
    for i in range(scale.shape[0]):
        if scale[i] <= 2:
            f_mel = np.append(f_mel,0.0255)
        if scale[i] >=3 and scale[i]<=5:
            f_mel = np.append(f_mel,0.155)
        if scale[i] == 6:
            f_mel = np.append(f_mel,0.305)
    return f_mel



def process_Binning(input_array,N):
    # Calculate the new shape of the array after reshaping
    new_shape = (input_array.shape[0],input_array.shape[1] // N, N)

    # Reshape the array to create chunks of size N
    reshaped_array = input_array[:,:input_array.shape[1] // N*N].reshape(new_shape)

    # Sum values along the second axis (axis=1)
    summed_array = np.sum(reshaped_array, axis=2)
    return summed_array

def get_minimum_sensitivity_threshold_normalization(binning,Intensity_measured,dr_phantom_simu):

    # dr_phantom_simu = process_Binning(dr_phantom_simu,binning)

    #Process binning on measured intensity
    binned_signal = process_Binning(Intensity_measured,binning)

    #Calculate SNR
    SNR = binned_signal.mean(axis=1)/binned_signal.std(axis=1)

    #Get coeff to scale simulations on experimental data
    Coeff = dr_phantom_simu.mean(axis=1)/binned_signal.mean(axis=1)

    #Noise level in measurement that corresponds to the minimum noise level
    # mup_min = binned_signal.std(axis=1)*Coeff
    mup_min = np.divide(binned_signal.mean(axis=1),SNR)*Coeff
    mup_min = np.ones(mup_min.shape)*np.max(mup_min)

    return mup_min, np.sqrt(mup_min), binned_signal.mean(axis=1)



def get_Photon_counts(dr,A,T_exp,wavelength_mm,S):
    return 50.34*A*T_exp*wavelength_mm*S*dr


## Load simulation data

path = main_path+ "simulations/output_lookup_table/"

#wavelength in nm
wavelength = 780

# Arrays thickness
skin_thickness_array = np.array([1,2,3,5])
# skin_thickness_array = np.array([1,2,3])

adipose_thickness_array = np.array([1,2,4,8,16])
muscle_thickness_array = np.array([3,7,10,13,20,25])

# Saturation array
SatO2_muscle_array = np.array([0.4,0.6,0.8])
SatO2_placenta_array = np.array([0.4,0.6,0.8])

#HbT
HbT_muscle_array = np.array([15,25,35,50])
HbT_placenta_array = np.array([15,25,35,50])

#Source detector separation simulation
SD_separations_mm_simulation = np.arange(10,90,10)



# Volume fraction of melanosome according to the color tones
f_melanosome = np.array([0.0255,0.155,0.305])

Tissue_sensitivity_indexes = np.zeros((skin_thickness_array.shape[0],
                                       adipose_thickness_array.shape[0],
                                       muscle_thickness_array.shape[0],
                                       SatO2_muscle_array.shape[0],
                                       SatO2_placenta_array.shape[0],
                                       f_melanosome.shape[0],
                                       HbT_muscle_array.shape[0],
                                       HbT_placenta_array.shape[0],
                                       SD_separations_mm_simulation.shape[0],
                                       4))

Diffuse_reflectance = np.zeros((skin_thickness_array.shape[0],
                                       adipose_thickness_array.shape[0],
                                       muscle_thickness_array.shape[0],
                                       SatO2_muscle_array.shape[0],
                                       SatO2_placenta_array.shape[0],
                                       f_melanosome.shape[0],
                                       HbT_muscle_array.shape[0],
                                       HbT_placenta_array.shape[0],
                                       SD_separations_mm_simulation.shape[0]))



#Data
for t_s in range(skin_thickness_array.shape[0]):
    for t_a in range(adipose_thickness_array.shape[0]):
        for t_m in range(muscle_thickness_array.shape[0]):
            for sat_m in range(SatO2_muscle_array.shape[0]):
                for sat_p in range(SatO2_placenta_array.shape[0]):
                    for f_mel in range(f_melanosome.shape[0]):
                        for HbT_m in range(HbT_muscle_array.shape[0]):
                            for HbT_p in range(HbT_placenta_array.shape[0]):

                                data = scipy.io.loadmat(path+
                                'out_St_muscle_' + str(SatO2_muscle_array[sat_m]) +
                                '_St_placenta_' + str(SatO2_placenta_array[sat_p]) +
                                '_Thick_skin_' +str(skin_thickness_array[t_s]) +
                                '_Thick_adipose_' +str(adipose_thickness_array[t_a]) +
                                '_Thick_muscle_' + str(muscle_thickness_array[t_m]) +
                                'f_mel' +str(f_melanosome[f_mel]) +
                                '_HbT_muscle_umol_'+str(HbT_muscle_array[HbT_m]) +
                                '_HbT_placenta_umol_'+str(HbT_placenta_array[HbT_p]) +'.mat')


                                Tissue_sensitivity_indexes[t_s,t_a,t_m,
                                sat_m,sat_p,f_mel,HbT_m,HbT_p,:,:] = data['Sensitivity_indexes']

                                Diffuse_reflectance[t_s,t_a,t_m,sat_m,sat_p,f_mel,HbT_m,HbT_p,:] = np.squeeze(data['Diffuse_reflectance'])




#Constants for converting diffuse reflectance in photons count
# Light source power (in W
Light_power_W = 1

# sensor area in mm2
Sensor_area_mm2 = 1e3

# wavelength in mm
wavelength_mm = wavelength*1e-3







# Load simulation of phantom measurements

# Load phantom
path = main_path + "simulations/"
data_phantom = scipy.io.loadmat(path+'out_Phantom.mat')
dr_phantom = data_phantom['Diffuse_reflectance']



#Load Monte Carlo noise of diffuse reflectance


files = glob.glob(main_path+"simulations/MC_noise_1GPU/*.mat")

dr_array = []
for f in files:
    data = scipy.io.loadmat(f)
    dr_array.append(data['Diffuse_reflectance'])

dr_array = np.squeeze(np.asarray(dr_array))
sigma_MC_noise_dr = dr_array.std(axis=0)






## Load experimental measurements with CYRIL and Mini CYRIL
path = main_path + "Phantoms/PhantomDataMiniCYRIL/"

SD_separations_cm_Mini_CYRIL = np.array([3,4,5])
SD_separations_mm_Mini_CYRIL = SD_separations_cm_Mini_CYRIL*10

integration_time_s_Mini_CYRIL = np.array([1,5,10])


#Size of the vector (integration time, SD separation, Time)
Intensity_measured_Mini_CYRIL = []

for j in range(integration_time_s_Mini_CYRIL.shape[0]):
    I_temp = []
    size = np.array([],dtype=int) #Size of the temporal vectors

    for i in range(SD_separations_cm_Mini_CYRIL.shape[0]):

            file = glob.glob(path+"*"+str(SD_separations_cm_Mini_CYRIL[i])+"cm_"+str(integration_time_s_Mini_CYRIL[j])+"s*")
            data = scipy.io.loadmat(file[0])

            Mini_CYRIL_wavelength = np.squeeze(data['Wavelengths'])

            Ref = np.squeeze(data['Ref'])
            Ref /= np.max(Ref)
            S = data['Spectra']

            # for t in range(S.shape[0]):
            #     S[i,:] /= Ref

            #ectract Ref!!!!!!!!!!!!!!!


            size = np.append(size,S.shape[0])
            I_temp.append(S)

    #crop the intensity to get the same vector lengths
    I_interp = np.zeros((SD_separations_cm_Mini_CYRIL.shape[0],size.min()))

    for i in range(SD_separations_cm_Mini_CYRIL.shape[0]):
        #Get intensity at specific wavelength
        for t in range(size.min()):
            I_interp[i,t] = interpolate.interp1d(Mini_CYRIL_wavelength,I_temp[i][t,:], kind='cubic')(wavelength)


    Intensity_measured_Mini_CYRIL.append(np.asarray(I_interp))
    print(I_interp.shape)


    #calculate SNR
    for binning in np.array([1,2,5,10]):
        binned_signal = process_Binning(I_interp,binning)

        #Get SNR_y
        SNR_y = binned_signal.mean(axis=1) / binned_signal.std(axis=1)


        print("ti ",integration_time_s_Mini_CYRIL[j],"binning", binning, "SNR ",SNR_y.mean(), "size signal: ",binned_signal.shape)


# plt.figure()
# for j in range(integration_time_s_Mini_CYRIL.shape[0]):
#
#     file = glob.glob(path+"*"+str(SD_separations_cm_Mini_CYRIL[0])+"cm_"+str(integration_time_s_Mini_CYRIL[j])+"s*")
#     data = scipy.io.loadmat(file[0])
#     plt.plot(Mini_CYRIL_wavelength,data['Spectra'].mean(axis=0),label=str(integration_time_s_Mini_CYRIL[j]))
# plt.legend(loc="best")
# plt.show()



# Load experimental measurements with CYRIL

#Path
path = main_path + "Phantoms/phantom_4mmL2L1L3__151223_113128/"

# Wavelength
CYRIL_wavelength = np.genfromtxt(main_path+"/raw_wavelengths.csv", delimiter=",")

#Integration time
integration_time_s_CYRIL = np.array([1])

#Source detector separation
detector_id = np.array([1,8,2,7,3,6,4,5]) #disposition of the detectors for SD_separations
SD_separations_mm_CYRIL = np.arange(10, detector_id.shape[0]*10+1, 10)


#Load intensity spectra (size: k, T, N), k: nb of detector, T: time, N: nb of wavelength
Intensity_measured_CYRIL = []

i = 0
for det_id in detector_id:
    file_path = glob.glob(path+"/Spectra/*"+str(det_id)+".csv")[0]
    I = np.genfromtxt(file_path, delimiter=",")


    #Get intensity at specific wavelength
    I_interp = np.array([])
    for t in range(I.shape[0]):
        I_interp = np.append(I_interp, interpolate.interp1d(CYRIL_wavelength,I[t,:], kind='cubic')(wavelength))

    Intensity_measured_CYRIL.append(I_interp)

#Convert list in ndarray
Intensity_measured_CYRIL = np.asarray(Intensity_measured_CYRIL)
temp = []
temp.append(Intensity_measured_CYRIL)
Intensity_measured_CYRIL = temp.copy()



## Save tissue sensitivity indexes


# # Use Mini CYRIL
# device_name = "Mini_CYRIL"
# integration_time_s = integration_time_s_Mini_CYRIL.copy()
# Intensity_measured = Intensity_measured_Mini_CYRIL.copy()
# SD_separations_mm = SD_separations_mm_Mini_CYRIL.copy()


#Use CYRIL
device_name = "CYRIL"
integration_time_s = integration_time_s_CYRIL.copy()
Intensity_measured = Intensity_measured_CYRIL.copy()
SD_separations_mm = SD_separations_mm_CYRIL.copy()

#Match detectors measured and simulated
id_det_simu = np.array([])
for i in SD_separations_mm:
    id_det_simu = np.append(id_det_simu,np.where(i==SD_separations_mm_simulation)[0].item())
id_det_simu = id_det_simu.astype(int)




#Save tissue sensitivity
np.savez(main_path+"simulations/tissue_sensitivity_"+device_name+".npz",
Tissue_sensitivity_indexes = Tissue_sensitivity_indexes[:,:,:,:,:,:,:,:,id_det_simu,:],
skin_thickness_array = skin_thickness_array,
adipose_thickness_array = adipose_thickness_array,
muscle_thickness_array = muscle_thickness_array,
SatO2_muscle_array = SatO2_muscle_array,
SatO2_placenta_array = SatO2_placenta_array,
f_melanosome = f_melanosome,
HbT_muscle_array = HbT_muscle_array,
HbT_placenta_array = HbT_placenta_array,
SD_separation_cm = SD_separations_mm/10
)


np.savez(main_path+"simulations/Diffuse_reflectance_"+device_name+".npz",
Diffuse_reflectance = Diffuse_reflectance[:,:,:,:,:,:,:,:,id_det_simu],
skin_thickness_array = skin_thickness_array,
adipose_thickness_array = adipose_thickness_array,
muscle_thickness_array = muscle_thickness_array,
SatO2_muscle_array = SatO2_muscle_array,
SatO2_placenta_array = SatO2_placenta_array,
f_melanosome = f_melanosome,
HbT_muscle_array = HbT_muscle_array,
HbT_placenta_array = HbT_placenta_array,
SD_separation_cm = SD_separations_mm/10
)

## Get minimum sensitivity threshold (normalize diffuse reflectance) - Mini CYRIL

lw = 5
ft = 18
ft_legend = 12

device = "Mini CYRIL"


plt.close('all')
fig = plt.figure(tight_layout=True)

gs = gridspec.GridSpec(1, 2)
ax1 = fig.add_subplot(gs[0,0])
ax3 = fig.add_subplot(gs[0,1])

# ax_simu = fig.add_subplot(gs[1,0])
# ax_mes = fig.add_subplot(gs[1,1])

id = 0
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

nb_binning = 10


label = []
SNR = []
mup_min_vec = []
mes_vec = []
simu_vec = []

if device == "CYRIL":
    integration_time_s = integration_time_s_CYRIL.copy()
    Intensity_measured = Intensity_measured_CYRIL.copy()
    for i in range(len(Intensity_measured)):
        Intensity_measured[i] = Intensity_measured[i][0:5,:]
    SD_separations_mm = SD_separations_mm_CYRIL.copy()
    SD_separations_mm = SD_separations_mm[0:5]

if device == "Mini CYRIL":
    integration_time_s = integration_time_s_Mini_CYRIL.copy()
    Intensity_measured = Intensity_measured_Mini_CYRIL.copy()
    SD_separations_mm = SD_separations_mm_Mini_CYRIL.copy()



#Match detectors measured and simulated
id_det_simu = np.array([])
for i in SD_separations_mm:
    id_det_simu = np.append(id_det_simu,np.where(i==SD_separations_mm_simulation)[0].item())
id_det_simu = id_det_simu.astype(int)



for i in range(integration_time_s.shape[0]):
    #Plot mu_p_phantom
    mup_min, sigma, binned_signal = get_minimum_sensitivity_threshold_normalization(1,Intensity_measured[i],dr_phantom[id_det_simu,:])

    mup_min_vec.append(mup_min)
    # mes_vec.append(binned_signal)
    # ax1.plot(SD_separations_mm,mup_min,color=colors[id],linewidth=lw,label=device+ " (IT = "+str(integration_time_s[i])+" s)")


    mup_min, sigma, binned_signal = get_minimum_sensitivity_threshold_normalization(nb_binning,Intensity_measured[i],dr_phantom[id_det_simu,:])
    mup_min_vec.append(mup_min)
    # mes_vec.append(binned_signal)
    # ax1.plot(SD_separations_mm,mup_min,color=colors[id+1],linewidth=lw,label=device+ " (IT = "+str(integration_time_s[i])+" s, temporal binning 10)")



    if device == "CYRIL":
        label.append("CYRIL (IT = "+str(integration_time_s[i])+" s)")
        SNR.append(Intensity_measured[i].mean(axis=1)/Intensity_measured[i].std(axis=1))
        mes_vec.append(Intensity_measured[i].std(axis=1))
        id +=1

    if device == "Mini CYRIL":
        SNR.append(Intensity_measured[i].mean(axis=1)/Intensity_measured[i].std(axis=1))
        mes_vec.append(Intensity_measured[i].std(axis=1))
        label.append("Mini CYRIL (IT = "+str(integration_time_s[i])+" s)")
        id +=1



    sig = process_Binning(Intensity_measured[i],nb_binning)

    if device == "CYRIL":
        label.append("CYRIL (IT = "+str(integration_time_s[i])+" s, temporal binning "+str(nb_binning)+")")
        SNR.append(sig.mean(axis=1)/sig.std(axis=1))
        mes_vec.append(sig.std(axis=1))
        id +=1

    if device == "Mini CYRIL":
        label.append("Mini CYRIL (IT = "+str(integration_time_s[i])+" s, temporal binning "+str(nb_binning)+")")
        SNR.append(sig.mean(axis=1)/sig.std(axis=1))
        mes_vec.append(sig.std(axis=1))
        id +=1






SNR = np.asarray(SNR).T
mup_min_vec = np.asarray(mup_min_vec).T
mes_vec = np.asarray(mes_vec).T

#Bar plot
width = 0.1
xbar = np.arange(SD_separations_mm.shape[0])


for i in range(SNR.shape[1]):
    offset = width * i

    # rects = ax1.bar(xbar + offset, mup_min_vec[:,i], width, label=label[i])
    rects = ax3.bar(xbar + offset, SNR[:,i], width, label=label[i])
    # rects = ax_mes.bar(xbar + offset, mes_vec[:,i], width, label=label[i])
    rects = ax1.bar(i,mup_min_vec[0,i], label=label[i])



    # ax_mes.plot(SD_separations_mm,mes_vec[:,i], label=label[i])


# ax_simu.plot(SD_separations_mm,dr_phantom[id_det_simu,:].mean(axis=1))

ax1.set_title("Minimum detectable diffuse reflectance of Mini CYRIL device",fontsize=ft)
ax1.legend(loc="best",fontsize=ft_legend)
# ax1.set_xlabel("Source detector separation (mm)",fontsize=ft)
ax1.set_ylabel("$\phi$ ($mm^{-2}.s^{-1}$)",fontsize=ft)
# ax1.set_xticks(np.arange(SD_separations_mm.shape[0])+5*width/2,np.array(["30","40","50"]),fontsize=ft)
ax1.set_xticks([],[])
ax1.tick_params(axis='both', which='major', labelsize=ft)
ax1.grid()

ax3.set_title("Signal to Noise Ratio (SNR) of Mini CYRIL device",fontsize=ft)
ax3.legend(loc="best",fontsize=ft_legend)
ax3.set_xlabel("Source detector separation (mm)",fontsize=ft)
ax3.set_ylabel("SNR",fontsize=ft)
ax3.set_xticks(np.arange(SD_separations_mm.shape[0])+5*width/2,np.array(["30","40","50"]),fontsize=ft)
ax3.grid()
ax3.tick_params(axis='both', which='major', labelsize=ft)

# ax_mes.set_yscale('log')
# ax_mes.set_xlabel("Source detector separation (mm)",fontsize=ft)
# ax_mes.set_ylabel("Source detector separation (mm)",fontsize=ft)
#
#
# ax_simu.set_yscale('log')
# ax_simu.set_xlabel("Source detector separation (mm)",fontsize=ft)

plt.show()


## Calculate the detection probability

# Use Mini CYRIL
device_name = "Mini CYRIL"
integration_time_s = integration_time_s_Mini_CYRIL.copy()
Intensity_measured = Intensity_measured_Mini_CYRIL.copy()
SD_separations_mm = SD_separations_mm_Mini_CYRIL.copy()


photon_count_model = 0

# #Use CYRIL
# device_name = "CYRIL"
# integration_time_s = integration_time_s_CYRIL.copy()
# Intensity_measured = Intensity_measured_CYRIL.copy()
# SD_separations_mm = SD_separations_mm_CYRIL.copy()


#Match detectors measured and simulated
id_det_simu = np.array([])
for i in SD_separations_mm:
    id_det_simu = np.append(id_det_simu,np.where(i==SD_separations_mm_simulation)[0].item())
id_det_simu = id_det_simu.astype(int)



for binning in np.array([1,10]):
    for id_int_time in range(integration_time_s.shape[0]):
        print("it",id_int_time,"binning",binning)

        #integration time
        int_time_s = integration_time_s[id_int_time]

        #Get mup min
        mu_p_min, sigma_p_phantom_detector,_ = get_minimum_sensitivity_threshold_normalization(binning,Intensity_measured[id_int_time],dr_phantom[id_det_simu,:])

        #log scale
        # mu_p_min = np.log10(mu_p_min)
        # sigma_p_phantom_detector = np.abs(np.log10(sigma_p_phantom_detector))


        #simulated data
        mu_p_detector = Diffuse_reflectance[:,:,:,:,:,:,:,:,id_det_simu]



        #Init output
        Detection_probability = np.zeros(mu_p_detector.shape)
        detector_signal = np.zeros(mu_p_detector.shape)
        MC_noise_proba = np.zeros(mu_p_detector.shape)


        for t_s in range(skin_thickness_array.shape[0]):
            for t_a in range(adipose_thickness_array.shape[0]):
                for t_m in range(muscle_thickness_array.shape[0]):
                    for sat_m in range(SatO2_muscle_array.shape[0]):
                        for sat_p in range(SatO2_placenta_array.shape[0]):
                            for f_mel in range(f_melanosome.shape[0]):
                                for HbT_m in range(HbT_muscle_array.shape[0]):
                                    for HbT_p in range(HbT_placenta_array.shape[0]):


                                        for d in range(SD_separations_mm.shape[0]):
                                            #Create normal distribution by addinf MC noise to diffuse
                                            #reflectance

                                            #Multiply by coeff and binning
                                            m_val = mu_p_detector[t_s,t_a,t_m,sat_m,sat_p,f_mel,HbT_m,HbT_p,d]*binning

                                            # #Do not multiply by binning (binning = mean)
                                            # m_val = mu_p_detector[t_s,t_a,t_m,sat_m,sat_p,f_mel,HbT_m,HbT_p,d]

                                            #log scale
                                            # m_val = np.log10(m_val)
                                            # if np.isnan(m_val) or np.isinf(m_val):
                                            #     continue





                                            # Detection probability
                                            sample = np.random.normal(loc = m_val,
                                                                    scale = sigma_p_phantom_detector[d],
                                                                    size = (100))
                                            # Process T-tests
                                            t,p = scipy.stats.ttest_1samp(sample,mu_p_min[d],
                                                                        alternative='greater')
                                            proba_detection = 1-p #Proba to accept alternative hypothesis





                                            #MC noise probability
                                            t,p = scipy.stats.ttest_1samp(sample,sigma_MC_noise_dr[d],
                                                                        alternative='greater')
                                            MC_proba = 1-p #Proba to accept alternative hypothesis







                                            # if m_val>mu_p_min[d]:
                                            #     proba_detection = 1
                                            # else:
                                            #     proba_detection = 0

                                            detector_signal[t_s,t_a,t_m,sat_m,sat_p,f_mel,HbT_m,HbT_p,d] = m_val
                                            Detection_probability[t_s,t_a,t_m,sat_m,sat_p,f_mel,HbT_m,HbT_p,d] = proba_detection
                                            MC_noise_proba[t_s,t_a,t_m,sat_m,sat_p,f_mel,HbT_m,HbT_p,d] = MC_proba



        #Save tissue sensitivity
        np.savez(main_path+"simulations/detection_proba/detection_probability_"+device_name+"_TI_"+str(int_time_s)+"s_binning_"+str(binning)+"_camera_model_"+str(photon_count_model)+".npz",
        Detection_probability = Detection_probability,
        MC_noise_proba = MC_noise_proba,
        detector_signal = detector_signal,
        Diffuse_reglectance = Diffuse_reflectance[:,:,:,:,:,:,:,:,id_det_simu],
        noise_level = sigma_p_phantom_detector,
        mu_p_min = mu_p_min,
        skin_thickness_array = skin_thickness_array,
        adipose_thickness_array = adipose_thickness_array,
        muscle_thickness_array = muscle_thickness_array,
        SatO2_muscle_array = SatO2_muscle_array,
        SatO2_placenta_array = SatO2_placenta_array,
        f_melanosome = f_melanosome,
        HbT_muscle_array = HbT_muscle_array,
        HbT_placenta_array = HbT_placenta_array,
        SD_separation_cm = SD_separations_mm/10
        )



##


sample = np.random.normal(loc = -10,
                        scale = 10,
                        size = (1000))
# Process T-tests
t,p = scipy.stats.ttest_1samp(sample,0,
                            alternative='greater')

print(1-p)
plt.close('all')
plt.hist(sample)
plt.show()