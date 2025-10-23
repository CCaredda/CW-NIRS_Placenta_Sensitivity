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
from numpy.linalg import lstsq

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



def read_thickness_data(main_path, data_patient,col=''):

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

    #remove PMS_ char in ID array
    for i in range(ID.shape[0]):
        ID[i] = ID[i][4:]
    ID = ID.astype(int)

    Fitzpatrick_scale = data['Fitzpatrick score'].values

    data_patient.Fitzpatrick_scale = np.zeros(data_patient.ID.shape)

    for i in range(ID.shape[0]):
        id = np.where(ID[i] == data_patient.ID)[0]
        if id.size == 0:
            continue
        data_patient.Fitzpatrick_scale[id] = Fitzpatrick_scale[i]



    return data_patient


def convert_in_Roman_Number(N):
    val = ""
    if N == 1:
        val = "I"
    if N == 2:
        val = "II"
    if N == 3:
        val = "III"
    if N == 4:
        val = "IV"
    if N == 5:
        val = "V"
    if N == 6:
        val = "VI"
    return val

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

    # dr_phantom_simu = (dr_phantom_simu,binning)

    #Process binning on measured intensity
    binned_signal = process_Binning(Intensity_measured,binning)

    #Calculate SNR
    SNR = binned_signal.mean(axis=1)/binned_signal.std(axis=1)

    #Get coeff to scale simulations on experimental data
    Coeff = np.divide(dr_phantom_simu,binned_signal.mean(axis=1))

    #Noise level in measurement that corresponds to the minimum noise level
    # mup_min = np.divide(binned_signal.mean(axis=1),SNR)*Coeff
    mup_min = binned_signal.std(axis=1)#/Coeff

    # mup_min = np.ones(mup_min.shape)*np.max(mup_min)

    return mup_min, mup_min, Coeff,
    # return mup_min, np.sqrt(mup_min), Coeff



def get_detection_proba(scaled_dr, sigma_p_phantom_detector, mu_p_min_detector):

    # Detection probability
    sample = np.random.normal(loc = scaled_dr,
                                scale = sigma_p_phantom_detector,
                                size = (1000))
    # Process T-tests
    t,p = scipy.stats.ttest_1samp(sample,mu_p_min_detector,
                                alternative='greater')
    return 1-p #Proba to accept alternative hypothesis

def SNR(signal):
    return np.mean(signal)/np.std(signal)

def SRS(Attenuation,SD_separations_in_mm,WAVELENGTHS,ext_coeffs_inv):

    #Calculate attenuation slope
    A_slope = (Attenuation[1, :] - Attenuation[0, :]) / (SD_separations_in_mm[1] - SD_separations_in_mm[0])


    #Calculate kmua
    h = 6.3*1e-4
    k_mua = (1 / (3 * (1 - (h * WAVELENGTHS))) )*  ((np.log(10) * A_slope -(2/np.mean(SD_separations_in_mm)))**2)

    #Get StO2
    kC = ext_coeffs_inv @ k_mua
    _StO2 = 100*(kC[1]/(kC[0]+kC[1]))

    if _StO2>100:
        _StO2 = 100
    if _StO2<0:
        _StO2 = 0
    return _StO2


## Read data measured by he clinical team

# Extract thickness from segmentation files
info_thickness = info_subject()
info_thickness = read_thickness_data(main_path, info_thickness)
info_thickness = read_thickness_data(main_path, info_thickness,col='.1')
info_thickness = read_thickness_data(main_path, info_thickness,col='.2')


ID_subjects = info_thickness.ID
Gestation_subjects = info_thickness.Gestation.astype(float)
t_skin = info_thickness.thickness_skin.astype(float)
t_muscle = info_thickness.thickness_muscle.astype(float)
t_adipose = info_thickness.thickness_adipose.astype(float)
Fitzpatrick = info_thickness.Fitzpatrick_scale.astype(float)


id_nan = ~ np.isnan(t_skin)
t_skin = t_skin[id_nan]
t_adipose = t_adipose[id_nan]
t_muscle = t_muscle[id_nan]
Fitzpatrick = Fitzpatrick[id_nan]
ID_subjects = ID_subjects[id_nan]
Gestation_subjects = Gestation_subjects[id_nan]

id_nan = ~ np.isnan(t_muscle)
t_skin = t_skin[id_nan]
t_adipose = t_adipose[id_nan]
t_muscle = t_muscle[id_nan]
Fitzpatrick = Fitzpatrick[id_nan]
ID_subjects = ID_subjects[id_nan]

Gestation_subjects = Gestation_subjects[id_nan]



dist_to_placenta = t_skin+t_adipose+t_muscle

# #remove subject with t_skin>0.5
# id_nan = t_skin<=0.3
# t_skin = t_skin[id_nan]
# t_adipose = t_adipose[id_nan]
# t_muscle = t_muscle[id_nan]
# Fitzpatrick = Fitzpatrick[id_nan]

#Convert Fitzpatrick scale in f_mel
f_mel = convert_Fitzpatrick_scale_f_mel(Fitzpatrick)




info = []
info.append(t_skin)
info.append(t_adipose)
info.append(t_muscle)
info.append(dist_to_placenta)
# info.append(Fitzpatrick)

print("fmel ",np.percentile(f_mel,25),np.percentile(f_mel,50),np.percentile(f_mel,75),np.percentile(f_mel,90))
print("skin (mm) ",int(np.percentile(t_skin,25)*10),int(np.percentile(t_skin,50)*10),int(np.percentile(t_skin,75)*10),int(np.percentile(t_skin,90)*10))
print("Adipose (mm)",int(np.percentile(t_adipose,25)*10),int(np.percentile(t_adipose,50)*10),int(np.percentile(t_adipose,75)*10),int(np.percentile(t_adipose,90)*10))
print("Muscle (mm)",int(np.percentile(t_muscle,25)*10),int(np.percentile(t_muscle,50)*10),int(np.percentile(t_muscle,75)*10),int(np.percentile(t_muscle,90)*10))

print("Distance to placenta (mm)",int(np.percentile(dist_to_placenta,25)*10),int(np.percentile(dist_to_placenta,50)*10),int(np.percentile(dist_to_placenta,75)*10),int(np.percentile(dist_to_placenta,90)*10))



#Save data
np.savetxt(main_path+"Subject_info/dist_to_placenta_cm.txt",dist_to_placenta)
np.savetxt(main_path+"Subject_info/t_skin_cm.txt",t_skin)
np.savetxt(main_path+"Subject_info/t_adipose_cm.txt",t_adipose)
np.savetxt(main_path+"Subject_info/t_muscle_cm.txt",t_muscle)

np.savetxt(main_path+"Subject_info/Fitzpatrick_scale.txt",Fitzpatrick)
np.savetxt(main_path+"Subject_info/Gestation_subjects_weeks.txt",Gestation_subjects)



## Plot data thickness data and gestational age

t_skin = np.loadtxt(main_path+"Subject_info/t_skin_cm.txt")
t_adipose = np.loadtxt(main_path+"Subject_info/t_adipose_cm.txt")
t_muscle = np.loadtxt(main_path+"Subject_info/t_muscle_cm.txt")

Fitzpatrick = np.loadtxt(main_path+"Subject_info/Fitzpatrick_scale.txt")
Gestation_subjects = np.loadtxt(main_path+"Subject_info/Gestation_subjects_weeks.txt")
dist_to_placenta = np.loadtxt(main_path+"Subject_info/dist_to_placenta_cm.txt")

info = []
info.append(t_skin)
info.append(t_adipose)
info.append(t_muscle)


ft = 23
plt.close('all')
plt.figure()
labels = np.array(["Skin","Adipose tissue", "Muscle"])
quantiles = [[0.25,0.70],[0.25,0.70],[0.25,0.70]]

plt.subplot(221)
plots = plt.violinplot(info,showmedians=True, showmeans=False)#,quantiles=quantiles)
# plots['cquantiles'].set_color('k')
# plots['cmeans'].set_color('r')
plots['cmedians'].set_color('r')
plt.xticks(np.arange(1,len(info)+1),labels,fontsize = ft)
plt.ylabel("Thickness (cm)",fontsize = ft)
plt.grid()
plt.title("Tissue thickness",fontsize = ft)
plt.yticks(fontsize=ft)


plt.subplot(222)
plt.title("Distance between skin and placenta",fontsize = ft)
plots = plt.violinplot(dist_to_placenta,showmedians=True, showmeans=False)#,quantiles=[0.25,0.75])
# plots['cquantiles'].set_color('k')
# plots['cmeans'].set_color('r')
plots['cmedians'].set_color('r')
plt.ylabel("Distance (cm)",fontsize = ft)
plt.xticks([],'')
plt.grid()
plt.yticks(fontsize=ft)


plt.subplot(223)
plt.title("Skin tones",fontsize = ft)
plots = plt.violinplot(Fitzpatrick,showmedians=True, showmeans=False)#,quantiles=[0.25,0.75])
# plots['cquantiles'].set_color('k')
# plots['cmeans'].set_color('r')
plots['cmedians'].set_color('r')
plt.ylabel("Fitzpatrick scale",fontsize = ft)
plt.xticks([],'')
plt.grid()
plt.yticks(fontsize=ft)

plt.subplot(224)
plt.title("Gestational age",fontsize = ft)
plots = plt.violinplot(Gestation_subjects,showmedians=True, showmeans=False)#,quantiles=[0.25,0.75])
# plots['cquantiles'].set_color('k')
# plots['cmeans'].set_color('r')
plots['cmedians'].set_color('r')
plt.ylabel("Gestational age (weeks)",fontsize = ft)
plt.xticks([],'')
plt.grid()
plt.yticks(fontsize=ft)
plt.show()



## Load experimental measurements with Mini CYRIL
path = main_path + "Phantoms/PhantomDataMiniCYRIL/"

SD_separations_cm = np.array([3,4,5])
SD_separations_mm = SD_separations_cm*10

integration_time_s = np.array([1,5,10])
wavelength = 780

binning_txt = np.array(["No binning", "binning 10"])

#Size of the vector (integration time, SD separation, Time)
Intensity_measured = []
Binned_signal = []
labels = []

for j in range(integration_time_s.shape[0]):
    I_temp = []
    size = np.array([],dtype=int) #Size of the temporal vectors

    for i in range(SD_separations_cm.shape[0]):

            file = glob.glob(path+"*"+str(SD_separations_cm[i])+"cm_"+str(integration_time_s[j])+"s*")
            data = scipy.io.loadmat(file[0])

            device_wavelength = np.squeeze(data['Wavelengths'])

            Ref = np.squeeze(data['Ref'])
            Ref /= np.max(Ref)
            S = data['Spectra']

            # for t in range(S.shape[0]):
            #     S[i,:] /= Ref

            #ectract Ref!!!!!!!!!!!!!!!


            size = np.append(size,S.shape[0])
            I_temp.append(S)

    #crop the intensity to get the same vector lengths
    I_interp = np.zeros((SD_separations_cm.shape[0],size.min()))

    for i in range(SD_separations_cm.shape[0]):
        #Get intensity at specific wavelength
        for t in range(size.min()):
            I_interp[i,t] = interpolate.interp1d(device_wavelength,I_temp[i][t,:], kind='cubic')(wavelength)


    Intensity_measured.append(np.asarray(I_interp))



    #calculate SNR
    for i_binning,binning in enumerate(np.array([1,10])):

        binned_signal = process_Binning(I_interp,binning)
        Binned_signal.append(binned_signal)
        labels.append("ti "+str(integration_time_s[j])+" s "+str(binning_txt[i_binning]))


ft=12
plt.rcParams.update({'font.size': ft})

plt.close('all')
plt.figure()
plt.subplot(231)
plt.title("Raw signal SD separation 3cm")
id = 0
plt.plot(Intensity_measured[0][id,:],'r',label=labels[0]+"\n SNR: "+str(int(SNR(Intensity_measured[0][id,:])))+
                                                "\n mean: "+str(int(np.mean(Intensity_measured[0][id,:]))) +
                                                "\n std: "+str(int(np.std(Intensity_measured[0][id,:]))))

plt.plot(Intensity_measured[2][id,:],'b',label=labels[4]+"\n SNR: "+str(int(SNR(Intensity_measured[2][id,:])))+
                                                "\n mean: "+str(int(np.mean(Intensity_measured[2][id,:]))) +
                                                "\n std: "+str(int(np.std(Intensity_measured[2][id,:]))))


plt.ylim(0,13000)
plt.legend(loc="best")
plt.ylabel("Intensity (a.u.)")
plt.xlabel("Temporal indexes")


plt.subplot(232)
plt.title("Raw signal SD separation 4cm")
id = 1
plt.plot(Intensity_measured[0][id,:],'r',label=labels[0]+"\n SNR: "+str(int(SNR(Intensity_measured[0][id,:])))+
                                                "\n mean: "+str(int(np.mean(Intensity_measured[0][id,:]))) +
                                                "\n std: "+str(int(np.std(Intensity_measured[0][id,:]))))

plt.plot(Intensity_measured[2][id,:],'b',label=labels[4]+"\n SNR: "+str(int(SNR(Intensity_measured[2][id,:])))+
                                                "\n mean: "+str(int(np.mean(Intensity_measured[2][id,:]))) +
                                                "\n std: "+str(int(np.std(Intensity_measured[2][id,:]))))

plt.legend(loc="best")
plt.ylim(0,13000)
plt.ylabel("Intensity (a.u.)")
plt.xlabel("Temporal indexes")


plt.subplot(233)
plt.title("Raw signal SD separation 5cm")
id = 2
plt.plot(Intensity_measured[0][id,:],'r',label=labels[0]+"\n SNR: "+str(int(SNR(Intensity_measured[0][id,:])))+
                                                "\n mean: "+str(int(np.mean(Intensity_measured[0][id,:]))) +
                                                "\n std: "+str(int(np.std(Intensity_measured[0][id,:]))))

plt.plot(Intensity_measured[2][id,:],'b',label=labels[4]+"\n SNR: "+str(int(SNR(Intensity_measured[2][id,:])))+
                                                "\n mean: "+str(int(np.mean(Intensity_measured[2][id,:]))) +
                                                "\n std: "+str(int(np.std(Intensity_measured[2][id,:]))))

plt.legend(loc="best")
plt.ylim(0,13000)
plt.ylabel("Intensity (a.u.)")
plt.xlabel("Temporal indexes")





plt.subplot(234)
plt.title("Binned signal SD separation 3cm")
id = 0
plt.plot(Binned_signal[1][id,:],'r',label=labels[1]+"\n SNR: "+str(int(SNR(Binned_signal[1][id,:])))+
                                                "\n mean: "+str(int(np.mean(Binned_signal[1][id,:]))) +
                                                "\n std: "+str(int(np.std(Binned_signal[1][id,:]))))

plt.plot(Binned_signal[5][id,:],'b',label=labels[5]+"\n SNR: "+str(int(SNR(Binned_signal[5][id,:])))+
                                                "\n mean: "+str(int(np.mean(Binned_signal[5][id,:]))) +
                                                "\n std: "+str(int(np.std(Binned_signal[5][id,:]))))
plt.ylim(0,130000)
plt.legend(loc="best")
plt.ylabel("Intensity (a.u.)")
plt.xlabel("Temporal indexes")

plt.subplot(235)
plt.title("Binned signal SD separation 4cm")
id = 1
plt.plot(Binned_signal[1][id,:],'r',label=labels[1]+"\n SNR: "+str(int(SNR(Binned_signal[1][id,:])))+
                                                "\n mean: "+str(int(np.mean(Binned_signal[1][id,:]))) +
                                                "\n std: "+str(int(np.std(Binned_signal[1][id,:]))))

plt.plot(Binned_signal[5][id,:],'b',label=labels[5]+"\n SNR: "+str(int(SNR(Binned_signal[5][id,:])))+
                                                "\n mean: "+str(int(np.mean(Binned_signal[5][id,:]))) +
                                                "\n std: "+str(int(np.std(Binned_signal[5][id,:]))))
plt.legend(loc="best")
plt.ylim(0,130000)
plt.ylabel("Intensity (a.u.)")
plt.xlabel("Temporal indexes")


plt.subplot(236)
plt.title("Binned signal SD separation 5cm")
id = 2
plt.plot(Binned_signal[1][id,:],'r',label=labels[1]+"\n SNR: "+str(int(SNR(Binned_signal[1][id,:])))+
                                                "\n mean: "+str(int(np.mean(Binned_signal[1][id,:]))) +
                                                "\n std: "+str(int(np.std(Binned_signal[1][id,:]))))

plt.plot(Binned_signal[5][id,:],'b',label=labels[5]+"\n SNR: "+str(int(SNR(Binned_signal[5][id,:])))+
                                                "\n mean: "+str(int(np.mean(Binned_signal[5][id,:]))) +
                                                "\n std: "+str(int(np.std(Binned_signal[5][id,:]))))
plt.legend(loc="best")
plt.ylim(0,130000)
plt.ylabel("Intensity (a.u.)")
plt.xlabel("Temporal indexes")

plt.show()






## Load and Plot tissue sensitivity data calculated with RedBird

wavelength = 890

#Load tissue sensitivity
path_data = main_path+"simulations/Redbird/data_article_"+str(wavelength)+"/"
HbT_placenta_array = np.array([15,25,35,50])
SD_separation_cm = np.array([3, 4, 5])

#23
ft = 23
ft_label = 23
plt.rcParams.update({'font.size': ft})

skin_thickness_subject_mm = np.array([1, 2, 3])
adipose_thickness_subject_mm = np.array([2, 4, 5])
muscle_thickness_subject_mm = np.array([7, 10, 12])
f_melanosome = np.array([2.55,15.5,30.5])/100
dist_to_plancenta_subjects_mm = skin_thickness_subject_mm + adipose_thickness_subject_mm + muscle_thickness_subject_mm


all_d = []

plt.close('all')
fig1 = plt.figure(tight_layout=True)
plt.suptitle("Placenta sensitivity at "+str(wavelength)+" nm",fontsize=ft)
# fig2 = plt.figure()

cmap = cm.plasma

# colors = ['k','w','w','w']
colors = ['k','k','k','k']

for id_mel in range(f_melanosome.shape[0]):
    for i in range(muscle_thickness_subject_mm.shape[0]):

        #Load placenta sensivitiy
        Placenta_sensitivity_redbird = np.zeros((HbT_placenta_array.shape[0], SD_separation_cm.shape[0]))
        for p,v_p in enumerate(HbT_placenta_array):

            data = scipy.io.loadmat(path_data+
            'out_St_muscle_0.6_St_placenta_0.8' +
            '_Thick_skin_' +str(skin_thickness_subject_mm[i]) +
            '_Thick_adipose_' +str(adipose_thickness_subject_mm[i]) +
            '_Thick_muscle_' + str(muscle_thickness_subject_mm[i]) +
            'f_mel' +str(f_melanosome[id_mel]) +
            '_HbT_muscle_umol_25'+
            '_HbT_placenta_umol_'+str(v_p) +'.mat')


            Placenta_sensitivity_redbird[p,:] = data['Sensitivity_indexes'][:,-1]


        #shape (fmel,detectors)

        all_d.append(Placenta_sensitivity_redbird)

        y = HbT_placenta_array
        x = SD_separation_cm.copy()
        T_x,T_y = np.meshgrid(x,y)

        ax1 = fig1.add_subplot(3,3,i+1+(3*(id_mel)))
        ax1.set_title("Placenta depth "+str(dist_to_plancenta_subjects_mm[i])+" mm",fontsize=ft)
        im = ax1.pcolor(T_x,T_y,Placenta_sensitivity_redbird*100, cmap=cmap, vmin=0, vmax=50)

        for x_id,xval in enumerate(x):
            for y_id,yval in enumerate(y):
                if(100*Placenta_sensitivity_redbird[y_id,x_id] < 10):
                    c="w"
                else:
                    c="k"

                ax1.text(xval-0.15,yval,str(int(10000*Placenta_sensitivity_redbird[y_id,x_id])/100),color=c,fontsize=ft)

        cb = fig1.colorbar(im, ax=ax1)
        cb.set_label("Placenta sensitivity (%)",fontsize=ft_label)
        ax1.set_xticks(x)  # Set x ticks to vector values
        ax1.set_yticks(y)  # Set y ticks to vector values
        ax1.set_xlabel("Source detector separation (cm)",fontsize=ft)
        ax1.set_ylabel("Placenta $C_{HbT}$ ($\mu$M)",fontsize=ft)

plt.show()

all_d = np.asarray(all_d)
print(wavelength, int(100*all_d.mean()))



## Influence of Skin in placenta sensitivity and detection proba

wavelength = 780

# Load simulation of phantom measurements
# Load phantom
path = main_path + "simulations/Redbird/"
data_phantom = scipy.io.loadmat(path+'Phantom_Data_'+str(wavelength)+'.mat')
dr_phantom = np.squeeze(data_phantom['DR_at_fiber_detector'])


#Load tissue sensitivity
path_data = main_path+"simulations/Redbird/data_article_subcutaneous/"
SD_separation_cm = np.array([3, 4, 5])

ft = 18
ft_label = 10
ft_text = 10
plt.rcParams.update({'font.size': ft_label})

skin_thickness_subject_mm = np.array([2, 2, 2, 2])
adipose_thickness_subject_mm = np.array([2, 4, 5, 7])
muscle_thickness_subject_mm = np.array([7, 10, 12, 17])
dist_to_plancenta_subjects_mm = skin_thickness_subject_mm + adipose_thickness_subject_mm + muscle_thickness_subject_mm

thickness_label = np.array([])
for i in range(skin_thickness_subject_mm.shape[0]):
    txt = "Placenta depth: "+str(dist_to_plancenta_subjects_mm[i])+" mm\n"+ "Skin: "+str(skin_thickness_subject_mm[i])+" mm\n"+"Adipose Tissue: "+str(adipose_thickness_subject_mm[i])+" mm\n"+"Muscle: "+str(muscle_thickness_subject_mm[i])+" mm"
    thickness_label = np.append(thickness_label,txt)

cmap = cm.plasma



binning = 1
int_time_s = 1
#Get mup min
mu_p_min, sigma_p_phantom_detector, Coeff = get_minimum_sensitivity_threshold_normalization(binning,Intensity_measured[np.where(int_time_s == integration_time_s)[0][0]],dr_phantom)





#Load Data
Placenta_sensitivity_redbird = np.zeros((skin_thickness_subject_mm.shape[0], SD_separation_cm.shape[0]))

Skin_sensitivity_redbird = np.zeros(Placenta_sensitivity_redbird.shape)
Adipose_sensitivity_redbird = np.zeros(Placenta_sensitivity_redbird.shape)
Muscle_sensitivity_redbird = np.zeros(Placenta_sensitivity_redbird.shape)
Detection_proba = np.zeros(Placenta_sensitivity_redbird.shape)

for i in range(skin_thickness_subject_mm.shape[0]):



    data = scipy.io.loadmat(path_data+
    'out_'+str(wavelength)+'_St_muscle_0.6_St_placenta_0.8' +
    '_Thick_skin_' +str(skin_thickness_subject_mm[i]) +
    '_Thick_adipose_' +str(adipose_thickness_subject_mm[i]) +
    '_Thick_muscle_' + str(muscle_thickness_subject_mm[i]) +
    'f_mel0.0255' +
    '_HbT_muscle_umol_35'+
    '_HbT_placenta_umol_35.mat')


    Skin_sensitivity_redbird[i,:] = data['Sensitivity_indexes'][:,0]
    Adipose_sensitivity_redbird[i,:] = data['Sensitivity_indexes'][:,1]

    Placenta_sensitivity_redbird[i,:] = data['Sensitivity_indexes'][:,3]
    Muscle_sensitivity_redbird[i,:] = data['Sensitivity_indexes'][:,2]

    val = np.squeeze(data['DR_at_fiber_detector'])/Coeff


    for d in range(SD_separation_cm.shape[0]):
        Detection_proba[i,d] = get_detection_proba(val[d], sigma_p_phantom_detector[d], mu_p_min[d])





skin_thickness_subject_mm = np.array([1, 2, 3])
adipose_thickness_subject_mm = np.array([5, 5, 5])
muscle_thickness_subject_mm = np.array([12, 12, 12])
dist_to_plancenta_subjects_mm2 = skin_thickness_subject_mm + adipose_thickness_subject_mm + muscle_thickness_subject_mm

thickness_label2 = np.array([])
for i in range(skin_thickness_subject_mm.shape[0]):
    txt = "Placenta depth: "+str(dist_to_plancenta_subjects_mm2[i])+" mm\n"+ "Skin: "+str(skin_thickness_subject_mm[i])+" mm\n"+"Adipose Tissue: "+str(adipose_thickness_subject_mm[i])+" mm\n"+"Muscle: "+str(muscle_thickness_subject_mm[i])+" mm"

    thickness_label2 = np.append(thickness_label2,txt)



#Load Data skin thickness influence
Placenta_sensitivity_redbird2 = np.zeros((skin_thickness_subject_mm.shape[0], SD_separation_cm.shape[0]))
Muscle_sensitivity_redbird2 = np.zeros(Placenta_sensitivity_redbird2.shape)
Detection_proba2 = np.zeros(Placenta_sensitivity_redbird2.shape)
Skin_sensitivity_redbird2 = np.zeros(Placenta_sensitivity_redbird2.shape)
Adipose_sensitivity_redbird2 = np.zeros(Placenta_sensitivity_redbird2.shape)

for i in range(skin_thickness_subject_mm.shape[0]):



    data = scipy.io.loadmat(path_data+
    'out_'+str(wavelength)+'_St_muscle_0.6_St_placenta_0.8' +
    '_Thick_skin_' +str(skin_thickness_subject_mm[i]) +
    '_Thick_adipose_' +str(adipose_thickness_subject_mm[i]) +
    '_Thick_muscle_' + str(muscle_thickness_subject_mm[i]) +
    'f_mel0.0255' +
    '_HbT_muscle_umol_35'+
    '_HbT_placenta_umol_35.mat')


    Skin_sensitivity_redbird2[i,:] = data['Sensitivity_indexes'][:,0]
    Adipose_sensitivity_redbird2[i,:] = data['Sensitivity_indexes'][:,1]

    Placenta_sensitivity_redbird2[i,:] = data['Sensitivity_indexes'][:,3]
    Muscle_sensitivity_redbird2[i,:] = data['Sensitivity_indexes'][:,2]

    val = np.squeeze(data['DR_at_fiber_detector'])/Coeff


    for d in range(SD_separation_cm.shape[0]):
        Detection_proba2[i,d] = get_detection_proba(val[d], sigma_p_phantom_detector[d], mu_p_min[d])



plt.close('all')
fig1 = plt.figure(tight_layout=True)
plt.suptitle("Effect of the skin thickness on detection probability",fontsize=ft)



#Plot 1
y = np.arange(dist_to_plancenta_subjects_mm.shape[0])
x = SD_separation_cm.copy()
T_x,T_y = np.meshgrid(x,y)
ax1 = fig1.add_subplot(1,2,1)
ax1.set_title("A - Detection probability at "+str(wavelength)+" nm for a fixed skin thickness",fontsize=ft)
im = ax1.pcolor(T_x,T_y,Detection_proba*100, cmap=cmap, vmin=0, vmax=100)

for x_id,xval in enumerate(x):
    for y_id,yval in enumerate(y):
        if Detection_proba[y_id,x_id]*100<40:
            colors = 'w'
        else:
            colors = 'k'
        ax1.text(xval-0.25,yval,str(int(10000*Detection_proba[y_id,x_id])/100),color=colors,fontsize=ft_text)

cb = fig1.colorbar(im, ax=ax1)
cb.set_label("Detection probability (%)",fontsize=ft_label)
ax1.set_xticks(x)  # Set x ticks to vector values
ax1.set_yticks(y)  # Set y ticks to vector values
ax1.set_yticklabels(thickness_label)
ax1.set_xlabel("Source detector separation (cm)",fontsize=ft)
# ax1.set_ylabel("Placenta depth (mm)",fontsize=ft)


#Plot 2
y = np.arange(dist_to_plancenta_subjects_mm2.shape[0])
x = SD_separation_cm.copy()
T_x,T_y = np.meshgrid(x,y)
ax1 = fig1.add_subplot(1,2,2)
ax1.set_title("A - Detection probability at "+str(wavelength)+" nm for varying skin thickness",fontsize=ft)
im = ax1.pcolor(T_x,T_y,Detection_proba2*100, cmap=cmap, vmin=0, vmax=100)

for x_id,xval in enumerate(x):
    for y_id,yval in enumerate(y):
        if Detection_proba2[y_id,x_id]*100<40:
            colors = 'w'
        else:
            colors = 'k'
        ax1.text(xval-0.25,yval,str(int(10000*Detection_proba2[y_id,x_id])/100),color=colors,fontsize=ft_text)

cb = fig1.colorbar(im, ax=ax1)
cb.set_label("Detection probability (%)",fontsize=ft_label)
ax1.set_xticks(x)  # Set x ticks to vector values
ax1.set_yticks(y)  # Set y ticks to vector values
ax1.set_yticklabels(thickness_label2)
ax1.set_xlabel("Source detector separation (cm)",fontsize=ft)
# ax1.set_ylabel("Placenta depth (mm)",fontsize=ft)
plt.show()




















fig1 = plt.figure(tight_layout=True)
plt.suptitle("Effect of the skin thickness on tissue sensitivity",fontsize=ft)


#Skin
#Plot 1
y = np.arange(dist_to_plancenta_subjects_mm.shape[0])
x = SD_separation_cm.copy()
T_x,T_y = np.meshgrid(x,y)
ax1 = fig1.add_subplot(2,4,1)
ax1.set_title("A - Skin sensitivity",fontsize=ft)
im = ax1.pcolor(T_x,T_y,Skin_sensitivity_redbird*100, cmap=cmap, vmin=0, vmax=100)

for x_id,xval in enumerate(x):
    for y_id,yval in enumerate(y):
        if Skin_sensitivity_redbird[y_id,x_id]*100<40:
            colors = 'w'
        else:
            colors = 'k'
        ax1.text(xval-0.25,yval,str(int(10000*Skin_sensitivity_redbird[y_id,x_id])/100),color=colors,fontsize=ft_text)

cb = fig1.colorbar(im, ax=ax1)
cb.set_label("Skin sensitivity (%)",fontsize=ft_label)
ax1.set_xticks(x)  # Set x ticks to vector values
ax1.set_yticks(y)  # Set y ticks to vector values
ax1.set_yticklabels(thickness_label)
ax1.set_xlabel("Source detector separation (cm)",fontsize=ft_label)
# ax1.set_ylabel("Placenta depth (mm)",fontsize=ft)


#Plot 2
y = np.arange(Placenta_sensitivity_redbird2.shape[0])
x = SD_separation_cm.copy()
T_x,T_y = np.meshgrid(x,y)
ax1 = fig1.add_subplot(2,4,2)
ax1.set_title("B - Skin sensitivity",fontsize=ft)
im = ax1.pcolor(T_x,T_y,Skin_sensitivity_redbird2*100, cmap=cmap, vmin=0, vmax=100)

for x_id,xval in enumerate(x):
    for y_id,yval in enumerate(y):
        if Skin_sensitivity_redbird2[y_id,x_id]*100<40:
            colors = 'w'
        else:
            colors = 'k'
        ax1.text(xval-0.25,yval,str(int(10000*Skin_sensitivity_redbird2[y_id,x_id])/100),color=colors,fontsize=ft_text)

cb = fig1.colorbar(im, ax=ax1)
cb.set_label("Skin sensitivity (%)",fontsize=ft_label)
ax1.set_xticks(x)  # Set x ticks to vector values
ax1.set_yticks(y)  # Set y ticks to vector values
ax1.set_yticklabels(thickness_label2)
ax1.set_xlabel("Source detector separation (cm)",fontsize=ft_label)
# ax1.set_ylabel("Placenta depth (mm)",fontsize=ft)







#Adipose
#Plot 1
y = np.arange(dist_to_plancenta_subjects_mm.shape[0])
x = SD_separation_cm.copy()
T_x,T_y = np.meshgrid(x,y)
ax1 = fig1.add_subplot(2,4,3)
ax1.set_title("C - Adipose sensitivity",fontsize=ft)
im = ax1.pcolor(T_x,T_y,Adipose_sensitivity_redbird*100, cmap=cmap, vmin=0, vmax=100)

for x_id,xval in enumerate(x):
    for y_id,yval in enumerate(y):
        if Adipose_sensitivity_redbird[y_id,x_id]*100<40:
            colors = 'w'
        else:
            colors = 'k'
        ax1.text(xval-0.25,yval,str(int(10000*Adipose_sensitivity_redbird[y_id,x_id])/100),color=colors,fontsize=ft_text)

cb = fig1.colorbar(im, ax=ax1)
cb.set_label("Adipose tissue sensitivity (%)",fontsize=ft_label)
ax1.set_xticks(x)  # Set x ticks to vector values
ax1.set_yticks(y)  # Set y ticks to vector values
ax1.set_yticklabels(thickness_label)
ax1.set_xlabel("Source detector separation (cm)",fontsize=ft_label)
# ax1.set_ylabel("Placenta depth (mm)",fontsize=ft)


#Plot 2
y = np.arange(Placenta_sensitivity_redbird2.shape[0])
x = SD_separation_cm.copy()
T_x,T_y = np.meshgrid(x,y)
ax1 = fig1.add_subplot(2,4,4)
ax1.set_title("D - Adipose sensitivity",fontsize=ft)
im = ax1.pcolor(T_x,T_y,Adipose_sensitivity_redbird2*100, cmap=cmap, vmin=0, vmax=100)

for x_id,xval in enumerate(x):
    for y_id,yval in enumerate(y):
        if Adipose_sensitivity_redbird2[y_id,x_id]*100<40:
            colors = 'w'
        else:
            colors = 'k'
        ax1.text(xval-0.25,yval,str(int(10000*Adipose_sensitivity_redbird2[y_id,x_id])/100),color=colors,fontsize=ft_text)

cb = fig1.colorbar(im, ax=ax1)
cb.set_label("Adispose tissue sensitivity (%)",fontsize=ft_label)
ax1.set_xticks(x)  # Set x ticks to vector values
ax1.set_yticks(y)  # Set y ticks to vector values
ax1.set_yticklabels(thickness_label2)
ax1.set_xlabel("Source detector separation (cm)",fontsize=ft_label)
# ax1.set_ylabel("Placenta depth (mm)",fontsize=ft)









#Muscle
#Plot 1
y = np.arange(dist_to_plancenta_subjects_mm.shape[0])
x = SD_separation_cm.copy()
T_x,T_y = np.meshgrid(x,y)
ax1 = fig1.add_subplot(2,4,5)
ax1.set_title("E - Muscle sensitivity",fontsize=ft)
im = ax1.pcolor(T_x,T_y,Muscle_sensitivity_redbird*100, cmap=cmap, vmin=0, vmax=100)

for x_id,xval in enumerate(x):
    for y_id,yval in enumerate(y):
        if Muscle_sensitivity_redbird[y_id,x_id]*100<40:
            colors = 'w'
        else:
            colors = 'k'
        ax1.text(xval-0.25,yval,str(int(10000*Muscle_sensitivity_redbird[y_id,x_id])/100),color=colors,fontsize=ft_text)

cb = fig1.colorbar(im, ax=ax1)
cb.set_label("Muscle sensitivity (%)",fontsize=ft_label)
ax1.set_xticks(x)  # Set x ticks to vector values
ax1.set_yticks(y)  # Set y ticks to vector values
ax1.set_yticklabels(thickness_label)
ax1.set_xlabel("Source detector separation (cm)",fontsize=ft_label)
# ax1.set_ylabel("Placenta depth (mm)",fontsize=ft)


#Plot 2
y = np.arange(Placenta_sensitivity_redbird2.shape[0])
x = SD_separation_cm.copy()
T_x,T_y = np.meshgrid(x,y)
ax1 = fig1.add_subplot(2,4,6)
ax1.set_title("F - Muscle sensitivity",fontsize=ft)
im = ax1.pcolor(T_x,T_y,Muscle_sensitivity_redbird2*100, cmap=cmap, vmin=0, vmax=100)

for x_id,xval in enumerate(x):
    for y_id,yval in enumerate(y):
        if Muscle_sensitivity_redbird2[y_id,x_id]*100<40:
            colors = 'w'
        else:
            colors = 'k'
        ax1.text(xval-0.25,yval,str(int(10000*Muscle_sensitivity_redbird2[y_id,x_id])/100),color=colors,fontsize=ft_text)

cb = fig1.colorbar(im, ax=ax1)
cb.set_label("Muscle sensitivity (%)",fontsize=ft_label)
ax1.set_xticks(x)  # Set x ticks to vector values
ax1.set_yticks(y)  # Set y ticks to vector values
ax1.set_yticklabels(thickness_label2)
ax1.set_xlabel("Source detector separation (cm)",fontsize=ft_label)
# ax1.set_ylabel("Placenta depth (mm)",fontsize=ft)












#Placenta
#Plot 1
y = np.arange(dist_to_plancenta_subjects_mm.shape[0])
x = SD_separation_cm.copy()
T_x,T_y = np.meshgrid(x,y)
ax1 = fig1.add_subplot(2,4,7)
ax1.set_title("G - Placenta sensitivity",fontsize=ft)
im = ax1.pcolor(T_x,T_y,Placenta_sensitivity_redbird*100, cmap=cmap, vmin=0, vmax=100)

for x_id,xval in enumerate(x):
    for y_id,yval in enumerate(y):
        if Placenta_sensitivity_redbird[y_id,x_id]*100<40:
            colors = 'w'
        else:
            colors = 'k'
        ax1.text(xval-0.25,yval,str(int(10000*Placenta_sensitivity_redbird[y_id,x_id])/100),color=colors,fontsize=ft_text)

cb = fig1.colorbar(im, ax=ax1)
cb.set_label("Placenta sensitivity (%)",fontsize=ft_label)
ax1.set_xticks(x)  # Set x ticks to vector values
ax1.set_yticks(y)  # Set y ticks to vector values
ax1.set_yticklabels(thickness_label)
ax1.set_xlabel("Source detector separation (cm)",fontsize=ft_label)
# ax1.set_ylabel("Placenta depth (mm)",fontsize=ft)


#Plot 2
y = np.arange(Placenta_sensitivity_redbird2.shape[0])
x = SD_separation_cm.copy()
T_x,T_y = np.meshgrid(x,y)
ax1 = fig1.add_subplot(2,4,8)
ax1.set_title("H - Placenta sensitivity",fontsize=ft)
im = ax1.pcolor(T_x,T_y,Placenta_sensitivity_redbird2*100, cmap=cmap, vmin=0, vmax=100)

for x_id,xval in enumerate(x):
    for y_id,yval in enumerate(y):
        if Placenta_sensitivity_redbird2[y_id,x_id]*100<40:
            colors = 'w'
        else:
            colors = 'k'
        ax1.text(xval-0.25,yval,str(int(10000*Placenta_sensitivity_redbird2[y_id,x_id])/100),color=colors,fontsize=ft_text)

cb = fig1.colorbar(im, ax=ax1)
cb.set_label("Placenta sensitivity (%)",fontsize=ft_label)
ax1.set_xticks(x)  # Set x ticks to vector values
ax1.set_yticks(y)  # Set y ticks to vector values
ax1.set_yticklabels(thickness_label2)
ax1.set_xlabel("Source detector separation (cm)",fontsize=ft_label)
# ax1.set_ylabel("Placenta depth (mm)",fontsize=ft)
plt.show()

## Calculate the detection probability


wavelength = 780
# Load simulation of phantom measurements
# Load phantom
path = main_path + "simulations/Redbird/"
data_phantom = scipy.io.loadmat(path+'Phantom_Data_'+str(wavelength)+'.mat')
dr_phantom = np.squeeze(data_phantom['DR_at_fiber_detector'])



ft = 18
ft_label = 18

plt.rcParams.update({'font.size': ft})


HbT_placenta_array = np.array([15,25,35,50])
skin_thickness_subject_mm = np.array([1, 2, 3])
adipose_thickness_subject_mm = np.array([2, 4, 5])
muscle_thickness_subject_mm = np.array([7, 10, 12])
f_melanosome = np.array([2.55,15.5,30.5])/100
dist_to_plancenta_subjects_mm = skin_thickness_subject_mm + adipose_thickness_subject_mm + muscle_thickness_subject_mm
SD_separation_cm = np.array([3, 4, 5])


colors = ['k','k','k','k']

#Set integration time and binning
int_time_s = 10
binning = 10



#Get mup min
mu_p_min, sigma_p_phantom_detector, Coeff = get_minimum_sensitivity_threshold_normalization(binning,Intensity_measured[np.where(int_time_s == integration_time_s)[0][0]],dr_phantom)



plt.close('all')
fig1 = plt.figure(tight_layout=True)
plt.suptitle("Detection probability at "+str(wavelength)+" nm - Integration time "+str(int_time_s)+" binning "+str(binning),fontsize=ft)

cmap = cm.plasma

# colors = ['k','w','w','w']
colors = ['k','k','k','k']


for id_mel in range(f_melanosome.shape[0]):
    for i in range(muscle_thickness_subject_mm.shape[0]):

        #Load placenta sensivitiy
        Detection_probability = np.zeros((HbT_placenta_array.shape[0], SD_separation_cm.shape[0]))


        for j,v_p in enumerate(HbT_placenta_array):

            data = scipy.io.loadmat(path+'data_article_'+str(wavelength) +
            '/out_St_muscle_0.6_St_placenta_0.8' +
            '_Thick_skin_' +str(skin_thickness_subject_mm[i]) +
            '_Thick_adipose_' +str(adipose_thickness_subject_mm[i]) +
            '_Thick_muscle_' + str(muscle_thickness_subject_mm[i]) +
            'f_mel' +str(f_melanosome[id_mel]) +
            '_HbT_muscle_umol_25'+
            '_HbT_placenta_umol_'+str(v_p) +'.mat')


            val = np.squeeze(data['DR_at_fiber_detector'])/Coeff


            for d in range(SD_separation_cm.shape[0]):
                Detection_probability[j,d] = get_detection_proba(val[d], sigma_p_phantom_detector[d], mu_p_min[d])



        y = HbT_placenta_array
        x = SD_separation_cm.copy()
        T_x,T_y = np.meshgrid(x,y)

        ax1 = fig1.add_subplot(3,3,i+1+(3*(id_mel)))
        ax1.set_title("Placenta depth "+str(dist_to_plancenta_subjects_mm[i])+" mm",fontsize=ft)
        im = ax1.pcolor(T_x,T_y,100*Detection_probability, cmap=cmap, vmin=0, vmax=100)

        for x_id,xval in enumerate(x):
            for y_id,yval in enumerate(y):
                ax1.text(xval-0.15,yval,str(int(10000*Detection_probability[y_id,x_id])/100),color=colors[i],fontsize=ft)

        cb = fig1.colorbar(im, ax=ax1)
        cb.set_label("Placenta sensitivity (%)",fontsize=ft_label)
        ax1.set_xticks(x)  # Set x ticks to vector values
        ax1.set_yticks(y)  # Set y ticks to vector values
        ax1.set_xlabel("Source detector separation (cm)",fontsize=ft)
        ax1.set_ylabel("Placenta $C_{HbT}$ ($\mu$M)",fontsize=ft)

plt.show()


## Plot detection proba, placenta sensitivity for all subjects



wavelength = 780
# Load simulation of phantom measurements
# Load phantom
path = main_path + "simulations/Redbird/"
data_phantom = scipy.io.loadmat(path+'Phantom_Data_'+str(wavelength)+'.mat')
dr_phantom = np.squeeze(data_phantom['DR_at_fiber_detector'])
SD_separation_cm = np.array([3, 4, 5])




ft = 23
ft_legend = 20
c = plt.rcParams['axes.prop_cycle'].by_key()['color']

plt.rcParams.update({'font.size': ft})

P_detection_m = []
P_scanning_m = []
x_ticks = []
x_labels = []
quantiles_vec = []
colors=[]
patches = []


data = scipy.io.loadmat(path+'Data_subjects_'+str(wavelength) + '.mat')
Diffuse_reflectance = data['DR_at_fiber_detector']
S_placenta = data['Sensitivity_indexes'][:,-1,:]
S_placenta_m = []
S_placenta_m.append(S_placenta[0,:])
S_placenta_m.append(S_placenta[1,:])
S_placenta_m.append(S_placenta[2,:])



for d in range(SD_separation_cm.shape[0]):
    for int_time_s in np.array([1,5,10]):
        for binning in np.array([1,10]):

            #Get mup min
            mu_p_min, sigma_p_phantom_detector, Coeff = get_minimum_sensitivity_threshold_normalization(binning,Intensity_measured[np.where(int_time_s == integration_time_s)[0][0]],dr_phantom)

            #init
            Detection_probability = np.zeros(S_placenta.shape[1],)
            Scanning_probability = np.zeros(S_placenta.shape[1],)


            for s in range(Diffuse_reflectance.shape[1]):
                val = Diffuse_reflectance[d,s]/Coeff[d]
                Detection_probability[s] = get_detection_proba(val, sigma_p_phantom_detector[d], mu_p_min[d])
                Scanning_probability[s] = Detection_probability[s]* S_placenta[d,s]




            # S_placenta_m.append(S_placenta[:,id_SatO2_muscle,id_SatO2_placenta,id_HbT_muscle,id_HbT_placenta, d])
            P_detection_m.append(Detection_probability)
            P_scanning_m.append(Scanning_probability)



for i in range(len(P_scanning_m)):
    quantiles_vec.append([0.25,0.75])

for i in range(3):
    for j in range(6):
        colors.append(c[j])

id = 0
for int_time_s in np.array([1,5,10]):
    for binning in np.array([1,10]):
        patches.append(mpatches.Patch(color=c[id], label="Integration time "+str(int_time_s)+" s, binning "+str(binning)))
        id += 1


plt.close('all')
fig = plt.figure(tight_layout=True)
gs = gridspec.GridSpec(4, 1, height_ratios=[1, 1, 0.1, 1])

ax = fig.add_subplot(gs[0,0])
ax.set_title("Placenta sensitivity",fontsize = ft)
plots = ax.violinplot(np.asarray(S_placenta_m).T*100, showmedians=False, showmeans=True,quantiles=[[0.25,0.75],[0.25,0.75],[0.25,0.75]])
plots['cmeans'].set_colors(c[id])
plots['cquantiles'].set_color('r')
ax.set_xticks(np.array([1,2,3]))
ax.set_xticklabels(np.array(["3","4","5"]),fontsize=ft)
ax.set_ylabel("Probability (%)",fontsize = ft)
ax.set_xlabel("Source-detector separation (cm)",fontsize = ft)

# for i in range(len(labels)-1):
#     ax.axvline(SD_cm.shape[0]+0.5+SD_cm.shape[0]*i)
ax.grid()
ax.set_ylim(0,100)


#Detection probability
ax = fig.add_subplot(gs[1,0])
ax.set_title("Detection probability",fontsize = ft)
plots = ax.violinplot(np.asarray(P_detection_m).T*100, showmedians=False, showmeans=True)#,quantiles=quantiles_vec)
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)
plots['cmeans'].set_colors(c[id])
# plots['cquantiles'].set_color('r')
ax.set_xticks(np.array([3.5,9.5,15.5]))
ax.set_xticklabels(np.array(["3","4","5"]),fontsize=ft)
ax.set_ylabel("Probability (%)",fontsize = ft)
ax.set_xlabel("Source-detector separation (cm)",fontsize = ft)
ax.axvline(6.5,0,100,color='k',linewidth=3)
ax.axvline(12.5,0,100,color='k',linewidth=3)
ax.grid()
ax.set_ylim(0,100)

# ax.legend(handles=patches,loc="upper right", bbox_to_anchor=(1.02, 1.0) ,fontsize = ft)


# Place one legend in the figure coordinates
# y = 0.37 is ~between 2nd and 3rd subplot → adjust as needed
fig.legend(handles=patches,
           loc="lower center", bbox_to_anchor=(0.5, 0.32),
           ncol=3, frameon=False,fontsize = ft_legend)




#Scanning probability
ax = fig.add_subplot(gs[3,0])
ax.set_title("Placenta scanning probability",fontsize = ft)
plots = ax.violinplot(np.asarray(P_scanning_m).T*100, showmedians=False, showmeans=True)#,quantiles=quantiles_vec)
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)
plots['cmeans'].set_colors(c[id])
# plots['cquantiles'].set_color('r')
ax.set_xticks(np.array([3.5,9.5,15.5]))
ax.set_xticklabels(np.array(["3","4","5"]),fontsize=ft)
ax.set_ylabel("Probability (%)",fontsize = ft)
ax.set_xlabel("Source-detector separation (cm)",fontsize = ft)
ax.axvline(6.5,0,100,color='k',linewidth=3)
ax.axvline(12.5,0,100,color='k',linewidth=3)
ax.grid()
ax.set_ylim(0,100)
# put those patched as legend-handles into the legend
# ax.legend(handles=patches,loc="best", borderaxespad=0. ,fontsize = ft)


        # id += 1
plt.show()

## Get minimum sensitivity threshold (normalize diffuse reflectance) - Mini CYRIL

wavelength = 780
# Load simulation of phantom measurements
# Load phantom
path = main_path + "simulations/Redbird/"
data_phantom = scipy.io.loadmat(path+'Phantom_Data_'+str(wavelength)+'.mat')
dr_phantom = np.squeeze(data_phantom['DR_at_fiber_detector'])



lw = 5
ft = 23
ft_legend = 20

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





for IT in np.array([1,10]):


    idx_IT = np.where(integration_time_s==IT)[0][0]
    #Plot mu_p_phantom
    mup_min, sigma, binned_signal = get_minimum_sensitivity_threshold_normalization(1,Intensity_measured[idx_IT],dr_phantom)

    mup_min_vec.append(mup_min)


    mup_min, sigma, binned_signal = get_minimum_sensitivity_threshold_normalization(nb_binning,Intensity_measured[idx_IT],dr_phantom)
    mup_min_vec.append(mup_min)


    SNR.append(Intensity_measured[idx_IT].mean(axis=1)/Intensity_measured[idx_IT].std(axis=1))
    mes_vec.append(Intensity_measured[idx_IT].std(axis=1))
    label.append("Mini CYRIL (IT = "+str(integration_time_s[idx_IT])+" s)")
    id +=1


    sig = process_Binning(Intensity_measured[idx_IT],nb_binning)
    label.append("Mini CYRIL (IT = "+str(integration_time_s[idx_IT])+" s, temporal binning "+str(nb_binning)+")")
    SNR.append(sig.mean(axis=1)/sig.std(axis=1))
    mes_vec.append(sig.std(axis=1))
    id +=1



SNR = np.asarray(SNR).T
mup_min_vec = np.asarray(mup_min_vec).T
mes_vec = np.asarray(mes_vec).T

#Bar plot
width = 0.2
xbar = np.arange(SD_separations_mm.shape[0])


for i in range(SNR.shape[1]):
    offset = width * i

    rects = ax1.bar(xbar + offset, mup_min_vec[:,i], width, label=label[i])
    rects = ax3.bar(xbar + offset, SNR[:,i], width, label=label[i])

    # rects = ax1.bar(xbar + offset, mes_vec[:,i], width, label=label[i])



    # ax_mes.plot(SD_separations_mm,mes_vec[:,i], label=label[i])



ax1.set_title("Minimum detectable diffuse reflectance of Mini CYRIL device",fontsize=ft)
ax1.legend(loc="best",fontsize=ft_legend)
ax1.set_xlabel("Source detector separation (mm)",fontsize=ft)
ax1.set_ylabel("$\phi_{min}$ (a. u.)",fontsize=ft)
ax1.set_xticks(np.arange(SD_separations_mm.shape[0])+3*width/2,np.array(["30","40","50"]),fontsize=ft)
ax1.tick_params(axis='both', which='major', labelsize=ft)
ax1.grid()

ax3.set_title("Signal to Noise Ratio (SNR) of Mini CYRIL device",fontsize=ft)
ax3.legend(loc="best",fontsize=ft_legend)
ax3.set_xlabel("Source detector separation (mm)",fontsize=ft)
ax3.set_ylabel("SNR",fontsize=ft)
ax3.set_xticks(np.arange(SD_separations_mm.shape[0])+3*width/2,np.array(["30","40","50"]),fontsize=ft)
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


## Study skin effect

path = main_path + "simulations/Redbird/"

#pos layers
Thickness_skin = 2
Thickness_adipose = 4
Thickness_muscle = 10

f_mel_array = np.array([0.0255,0.305])
SD_separation_cm_array = np.array([3, 4, 5])
SD_separation_id = 0


data = scipy.io.loadmat(path+"Sensitivity_proba_fmel_"+str(f_mel_array[0])+".mat")
srcpos = np.squeeze(data['srcpos'])
detpos = np.squeeze(data['detpos'])
thickness_layers_mm = np.squeeze(data['thickness_layers_mm'])

depth = 30

# levels = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100]
levels = [1e-3,5e-3,1e-2,5e-2,1e-1,5e-1,1,5,10,50,100]


ft = 23
ft_txt = 20
lw = 3
plt.close('all')
plt.figure()

plt.rcParams.update({'font.size': ft})


vec_map_interp = []

for f_mel in range(f_mel_array.shape[0]):

    #Load data
    data = scipy.io.loadmat(path+"Sensitivity_proba_fmel_"+str(f_mel_array[f_mel])+".mat")
    sensitivity_map = data['Sensitivity_proba']
    #Select sensitivity for the correct SD separation
    sensitivity_map = sensitivity_map[:,:,:,SD_separation_id]
    #Normalize to get proba densitity
    sensitivity_map = sensitivity_map/sensitivity_map.max()


    #Convert in percentage
    sensitivity_map *= 100

    # Source pos
    volume_square_size = 200
    src_pos = np.array([srcpos[0] ,srcpos[1], 0])

    #Define detector
    det_pos = np.array([detpos[SD_separation_id,0], detpos[SD_separation_id,1], 0])



    # Create cross section
    # map = sensitivity_map[det_pos[0],src_pos[1]-10:det_pos[1]+10,0:22].copy()
    map = sensitivity_map[src_pos[1]-10:det_pos[1]+10,det_pos[0],0:depth].copy()
    map = map.T

    #interpolate
    reso = 1
    x = np.arange(0,map.shape[1]*reso,reso)
    y = np.arange(0,map.shape[0]*reso,reso)
    interp_map = interpolate.RegularGridInterpolator((y, x),map,bounds_error=False, fill_value=None)

    x_interp = np.linspace(x[0],x[-1],x.shape[0]*4)
    y_interp = np.linspace(y[0],y[-1],y.shape[0]*4)

    map_interp = np.zeros((y_interp.shape[0],x_interp.shape[0]))
    for i in range(x_interp.shape[0]):
        for j in range(y_interp.shape[0]):
            map_interp[j,i] = interp_map((y_interp[j], x_interp[i]))

    vec_map_interp.append(map_interp)
    p = plt.subplot(1,2,f_mel+1)
    plt.title("Melanosome volume fraction "+str(f_mel_array[f_mel]*100)+"%",fontsize=ft)
    im = plt.contourf(x_interp, y_interp,map_interp,levels = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100], locator=ticker.LogLocator(),cmap='plasma')

    plt.plot(x_interp,thickness_layers_mm[0]*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)
    plt.plot(x_interp,thickness_layers_mm[1]*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)
    plt.plot(x_interp,thickness_layers_mm[2]*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)

    plt.plot(10,0,'ko',markersize=12,label="Source")
    plt.plot(10+SD_separation_cm_array[SD_separation_id]*10,0,'go',markersize=12,label="Detector")

    plt.text(1,thickness_layers_mm[0]-0.2,"Skin",fontsize=ft_txt)
    plt.text(1,thickness_layers_mm[1]-0.2,"Adipose tissue",fontsize=ft_txt)
    plt.text(1,thickness_layers_mm[2]-0.2,"Muscle",fontsize=ft_txt)
    plt.text(1,thickness_layers_mm[2]+0.5,"Placenta",fontsize=ft_txt)





    c = plt.colorbar(im)
    c.set_label("Sensitivity probability (%)",fontsize=ft)
    p.invert_yaxis()
    # plt.xticks(x)
    # plt.yticks(y)
    plt.xlabel("Tissue width (mm)",fontsize=ft)
    plt.ylabel("Tissue depth (mm)",fontsize=ft)
    plt.legend(loc="best",fontsize=ft)


plt.show()



## Calculate and plot SRS

# Load phantom
path = main_path + "simulations/Redbird/"
data_phantom = scipy.io.loadmat(path+'Phantom_Data_780.mat')
dr_phantom = np.squeeze(data_phantom['DR_at_fiber_detector'])


#Set integration time and binning
int_time_s = 10
binning = 10



#Get mup min
mu_p_min, sigma_p_phantom_detector, Coeff = get_minimum_sensitivity_threshold_normalization(binning,Intensity_measured[np.where(int_time_s == integration_time_s)[0][0]],dr_phantom)


#Extinction coeff shape (wavelengths,nb_chrom)
# chromophores are Hb, HbO2
extinction_coefficients = np.array([ [0.000110593822843823, 7.54002331002331e-05],
                                        [7.98426573426574e-05, 9.14371212121212e-05],
                                        [7.77872960372961e-05, 0.000100791258741259],
                                        [7.77082750582751e-05, 0.000105600582750583],
                                        [7.81447552447553e-05, 0.000109727972027972],
                                        [8.63153846153846e-05, 0.000122350815850816]])

# Get MBLL Matrix
ext_coeffs_inv = np.linalg.pinv(extinction_coefficients)

#Wavelengths
WAVELENGTHS = np.array([780, 810, 830, 840, 850, 890])

#Source detector separation
SD_separations_in_mm = np.array([30,40,50])
id_comb_34 = np.array([0,1]) #3, 4 cm
id_comb_45 = np.array([1,2]) #4, 5 cm

#SatO2
SatO2_muscle = np.array([0.4, 0.4, 0.5, 0.7, 0.8, 0.8])
SatO2_placenta = np.array([0.4, 0.8, 0.7, 0.5, 0.4, 0.8])




#fmelanosome
f_melanosome = np.array([2.55,15.5,30.5])/100

#Subject tissue thickness
skin_thickness_subject_mm = np.array([1, 2, 3])
adipose_thickness_subject_mm = np.array([2, 4, 5])
muscle_thickness_subject_mm = np.array([7, 10, 12])
dist_to_plancenta_subjects_mm = skin_thickness_subject_mm + adipose_thickness_subject_mm + muscle_thickness_subject_mm

#Init output
SatO2_array = np.zeros((3,SatO2_muscle.shape[0],f_melanosome.shape[0]))
detection_proba = np.zeros(SatO2_array.shape)


for id_subject in range(1,4):
    for S in range(SatO2_muscle.shape[0]):
        for f_mel in range(f_melanosome.shape[0]):
            #Load mat file
            data = scipy.io.loadmat(main_path + "simulations/Redbird/StO2_LUT/out_subject_" +
                                    str(id_subject)+"St_muscle_"+str(SatO2_muscle[S]) +
                                    "_St_placenta_"+str(SatO2_placenta[S])+ "_fmel_" +
                                    str(f_melanosome[f_mel])+".mat")

            #Load intensity
            I = np.squeeze(data['DR_at_fiber_detector'])

            #Select detector comb
            I = I[id_comb_34,:]

            #Calculate detection proba at 780 nm for detector at 3 cm
            det1 = get_detection_proba(I[id_comb_34[0],0]/Coeff[id_comb_34[0]], sigma_p_phantom_detector[id_comb_34[0]], mu_p_min[id_comb_34[0]])
            det2 = get_detection_proba(I[id_comb_34[1],0]/Coeff[id_comb_34[1]], sigma_p_phantom_detector[id_comb_34[1]], mu_p_min[id_comb_34[1]])

            if det1<0.99 or det2<0.99:
                detection_proba[id_subject-1, S, f_mel] = 0
            else:
                detection_proba[id_subject-1, S, f_mel] = 1

            #Compute Attenuation
            Attenuation = np.log10(1/I)

            #Process SRS
            _StO2 = SRS(Attenuation[id_comb_34,:],
                        SD_separations_in_mm[id_comb_34],WAVELENGTHS,ext_coeffs_inv)

            SatO2_array[id_subject-1, S, f_mel] = _StO2




#Set to nan if detection proba is below 99%
SatO2_array[detection_proba<0.99] = np.nan



#SatO2 homogeneous tissue
SatO2_homo = np.array([0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
SatO2_array_homo = np.zeros((2,SatO2_homo.shape[0]))
detection_proba_homo = np.zeros(SatO2_array_homo.shape)
id_comb_34_homo = np.array([0,2])
id_comb_45_homo = np.array([2,4])


for S in range(SatO2_homo.shape[0]):
    #Load mat file
    data = scipy.io.loadmat(main_path + "simulations/Redbird/StO2_semi_ininite/St_" +str(SatO2_homo[S])+".mat")

    #Load intensity
    I = np.squeeze(data['DR_at_fiber_detector'])

    #Compute Attenuation
    Attenuation = np.log10(1/I)


    #Calculate detection proba at 780 nm for detector at 3 cm
    detection_proba_homo[0,S] = get_detection_proba(I[id_comb_34_homo[0],0]/Coeff[id_comb_34_homo[0]], sigma_p_phantom_detector[id_comb_34_homo[0]], mu_p_min[id_comb_34_homo[0]])
    detection_proba_homo[1,S] = get_detection_proba(I[id_comb_34_homo[1],0]/Coeff[id_comb_34_homo[1]], sigma_p_phantom_detector[id_comb_34_homo[1]], mu_p_min[id_comb_34_homo[1]])

    #Process SRS
    _StO2_SD_34 = SRS(Attenuation[id_comb_34_homo,:], np.array([30,40]),
                        WAVELENGTHS,ext_coeffs_inv)
    _StO2_SD_45 = SRS(Attenuation[id_comb_45_homo,:], np.array([40,50]),
                        WAVELENGTHS,ext_coeffs_inv)

    SatO2_array_homo[0,S] = _StO2_SD_34
    SatO2_array_homo[1,S] = _StO2_SD_45



#Set to nan if detection proba is below 95%
SatO2_array_homo[detection_proba_homo<0.99] = np.nan


# plot
ft = 20
ft_label = ft
plt.rcParams.update({'font.size': ft})

val_thresh = 50


SatO2_homo_display = SatO2_array_homo[0,:]
SatO2_homo_display = SatO2_homo_display.reshape(SatO2_homo_display.shape[0], 1)

SatO2_array_display = SatO2_array.copy()

plt.close('all')
fig1, axes = plt.subplots(4, 1, constrained_layout=True)
SD = (SD_separations_in_mm[id_comb_34]/10).astype(int)

SD_sep = "(SD separation: "+str(SD[0])+"-"+str(SD[1])+" cm)"


# ----- Display SatO2 for semi-infinite slab -----
axes[0].set_title("$PltO_2$ estimation for a semi-infinite slab "+SD_sep, fontsize=ft)

# Fix shape for pcolor
x = np.arange(SatO2_homo_display.shape[1] + 1)  # Columns +1
y = np.arange(SatO2_homo_display.shape[0] + 1)  # Rows +1

im = axes[0].pcolor(y, x, SatO2_homo_display.T, vmin=0, vmax=100, cmap=cmap)

# Custom x labels
x_positions = np.arange(SatO2_homo_display.shape[0]) + 0.5
x_labels = []
for s in SatO2_homo:
    x_labels.append(str(int(s*100))+"%")

axes[0].set_yticks([])  # Remove y-ticks
axes[0].set_xticks(x_positions)  # Set positions for labels
axes[0].set_xticklabels(x_labels, fontsize=ft)
# axes[0].set_xlabel("$SatO_2$ value in the semi-infinite slab (%)", fontsize=ft)

# Add text annotations
for x_id, xval in enumerate(x_positions):
    val = SatO2_homo_display[x_id, 0]
    if not np.isnan(val):
        if val<val_thresh:
            col = 'w'
        else:
            col = 'k'
        axes[0].text(xval-0.05, 0.5, str(int(val)), color=col, fontsize=ft)




#xlabel
x_labels = np.array([])
for S in range(SatO2_muscle.shape[0]):
    x_labels = np.append(x_labels, "Muscle: "+str(int(SatO2_muscle[S]*100))+"%"+"\n"+
                                    "Placenta: "+str(int(SatO2_placenta[S]*100))+"%")

#ylabel
y_labels = np.array([])
for S in range(f_melanosome.shape[0]):
    y_labels = np.append(y_labels, str((f_melanosome[S]*100))+"%")


for subject_id in range(SatO2_array_display.shape[0]):
    ax = axes[subject_id + 1]  # Access subplot

    # Define meshgrid correctly
    x = np.arange(x_labels.shape[0] + 1)
    y = np.arange(y_labels.shape[0] + 1)
    T_x, T_y = np.meshgrid(x, y)

    im = ax.pcolor(T_x, T_y, SatO2_array_display[subject_id, :, :].T, vmin=0, vmax=100, cmap=cmap)

    ax.set_xticks(np.arange(len(x_labels)) + 0.5)
    ax.set_xticklabels(x_labels, fontsize=ft)

    ax.set_yticks(np.arange(len(y_labels)) + 0.5)
    ax.set_yticklabels(y_labels, fontsize=ft)

    ax.set_title("$PltO_2$ estimation for a placenta depth of "+str(dist_to_plancenta_subjects_mm[subject_id]) +" mm "+SD_sep, fontsize=ft)

    # Add text annotations
    for x_id in range(len(x_labels)):
        for y_id in range(len(y_labels)):
            val = SatO2_array_display[subject_id, x_id, y_id]
            if not np.isnan(val):
                if val<val_thresh:
                    col = 'w'
                else:
                    col = 'k'
                ax.text(x_id + 0.5, y_id +0.4, str(int(val)), color=col, fontsize=ft)

# Add shared labels
fig1.supylabel("Melanosome volume fraction (%)", fontsize=ft)
fig1.supxlabel("Placental tissue oxygenation (%)", fontsize=ft)

# Fix colorbar expansion issue
cbar = fig1.colorbar(im, ax=axes, orientation='vertical', fraction=0.02, pad=0.04, aspect=30)
cbar.set_label("$PltO_2$ (%)", fontsize=ft)

plt.show()







ms = 15

# Plot PltO2 for semi infinite slab

plt.figure()
plt.title("Placental tissue oxygenation estimation in semi-infinite slab with SRS algorithm",fontsize=ft)
plt.plot(SatO2_homo*100,SatO2_homo*100,"k",label="Ground truth")

plt.plot(SatO2_homo*100,SatO2_array_homo[0,:],"ko",label="Source-detector separation 3-4 cm",markersize=ms)
plt.plot(SatO2_homo*100,SatO2_array_homo[1,:],"r*",label="Source-detector separation 4-5 cm",markersize=ms)
plt.xlabel("Expected $PltO_2$ (%)",fontsize=ft)
plt.ylabel("Estimated $PltO_2$ (%)",fontsize=ft)
# plt.xlim(0,100)
# plt.ylim(0,100)

plt.legend(loc="best",fontsize=ft)
plt.grid()

plt.show()




#Errors
e_homo = (SatO2_array_homo[0,:] - SatO2_homo*100)


##

SD_separations_in_mm = np.array([30,35,40,45,50])
SD_separations_in_cm = np.array(["3", "3.5","4", "4.5", "5"])

id_comb = np.array([[0,1],[3,4],[0,2],[2,4]])
colors = np.array(["ko","ro","go","mo"])

#SatO2 homogeneous tissue
SatO2_homo = np.array([0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
SatO2_array_homo = np.zeros((id_comb.shape[0],SatO2_homo.shape[0]))


for S in range(SatO2_homo.shape[0]):
    #Load mat file
    data = scipy.io.loadmat(main_path + "simulations/Redbird/StO2_semi_ininite/St_" +str(SatO2_homo[S])+".mat")

    #Load intensity
    I = np.squeeze(data['DR_at_fiber_detector'])

    #Compute Attenuation
    Attenuation = np.log10(1/I)


    #Process SRS
    for id in range(id_comb.shape[0]):
        _StO2 = SRS(Attenuation[id_comb[id,:],:], SD_separations_in_mm[id_comb[id,:]],
                        WAVELENGTHS,ext_coeffs_inv)
        SatO2_array_homo[id,S] = _StO2


plt.figure()
plt.title("Placental tissue oxygenation estimation in semi-infinite slab with SRS algorithm",fontsize=ft)
plt.plot(SatO2_homo*100,SatO2_homo*100,"k",label="Ground truth")

for id in range(id_comb.shape[0]):
    plt.plot(SatO2_homo*100,SatO2_array_homo[id,:],colors[id],label="Source-detector separation "+SD_separations_in_cm[id_comb[id,0]]+"-"+SD_separations_in_cm[id_comb[id,1]]+" cm",markersize=ms)

plt.xlabel("Expected $PltO_2$ (%)",fontsize=ft)
plt.ylabel("Estimated $PltO_2$ (%)",fontsize=ft)
# plt.xlim(0,100)
# plt.ylim(0,100)

plt.legend(loc="best",fontsize=ft)
plt.grid()

plt.show()


## Compare MCX and Redbird values


path = main_path + "simulations/Redbird/compare_MCX_Redbird/"


thickness_layers_mm = np.array([3, 5, 12])

SD_separation_cm_array = np.array([3, 4, 5])
SD_separation_id = 0


data_MCX = scipy.io.loadmat(path+"MCX.mat")
MCX_DR = data_MCX['Diffuse_reflectance']
MCX_Sensitivity_profile = data_MCX['Sensitivity_profile']
MCX_Sensitivity_index = data_MCX['Sensitivity_indexes']


data_Redbird = scipy.io.loadmat(path+"Redbird.mat")
Redbird_DR = data_Redbird['DR_at_fiber_detector']
Redbird_Sensitivity_profile = data_Redbird['sensitivity_profile']
Redbird_Sensitivity_index = data_Redbird['Sensitivity_indexes']

src_pos = np.squeeze(data_Redbird['src_pos'])
det_pos = np.squeeze(data_Redbird['det_pos'])
det_pos = det_pos[SD_separation_id]


depth = 30
S = MCX_Sensitivity_profile[:,:,:,SD_separation_id]
S = S/S.max()
map_MCX = S[det_pos[0],src_pos[1]-10:det_pos[1]+10,0:depth].copy()
map_MCX = 100*map_MCX.T
map_MCX= map_MCX[1:,:] # remove first row for redbird compa (first row is air)

S = Redbird_Sensitivity_profile[:,:,:,SD_separation_id]
S = S/S.max()
map_Redbird = S[src_pos[1]-10:det_pos[1]+10,det_pos[0],0:depth].copy()
map_Redbird = 100*map_Redbird.T
map_Redbird  = map_Redbird[0:-1,:]

# levels = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100]
levels = [1e-3,5e-3,1e-2,5e-2,1e-1,5e-1,1,5,10,50,100]



ft = 23
ft_txt = 20
lw = 3
plt.close('all')
plt.figure()

plt.rcParams.update({'font.size': ft})




#interpolate
reso = 1
x = np.arange(0,map_Redbird.shape[1]*reso,reso)
y = np.arange(0,map_Redbird.shape[0]*reso,reso)
interp_map_Redbird = interpolate.RegularGridInterpolator((y, x),map_Redbird,bounds_error=False, fill_value=None)
interp_map_MCX = interpolate.RegularGridInterpolator((y, x),map_MCX,bounds_error=False, fill_value=None)

x_interp = np.linspace(x[0],x[-1],x.shape[0]*4)
y_interp = np.linspace(y[0],y[-1],y.shape[0]*4)


map_interp_Redbird = np.zeros((y_interp.shape[0],x_interp.shape[0]))
map_interp_MCX = np.zeros((y_interp.shape[0],x_interp.shape[0]))

for i in range(x_interp.shape[0]):
    for j in range(y_interp.shape[0]):
        map_interp_Redbird[j,i] = interp_map_Redbird((y_interp[j], x_interp[i]))
        map_interp_MCX[j,i] = interp_map_MCX((y_interp[j], x_interp[i]))



vec_map = [map_interp_Redbird, map_interp_MCX, map_interp_Redbird - map_interp_MCX]
title = ["A - Redbird", "B - MCX", "C - Difference between\nRedbird and MCX"]

plt.close('all')
plt.figure()

levels_contours = [1e-5,1e-4,1e-3,1e-2,1e-1,1,30,50,100]

for i in range(2):
    p=plt.subplot(1,3,i+1)
    plt.title(title[i],fontsize=ft)
    im = plt.contourf(x_interp, y_interp,vec_map[i],levels = levels_contours, locator=ticker.LogLocator(),cmap='plasma')

    plt.plot(x_interp,thickness_layers_mm[0]*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)
    plt.plot(x_interp,thickness_layers_mm[1]*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)
    plt.plot(x_interp,thickness_layers_mm[2]*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)

    plt.plot(10,0,'ko',markersize=12,label="Source")
    plt.plot(10+SD_separation_cm_array[SD_separation_id]*10,0,'go',markersize=12,label="Detector")

    plt.text(1,thickness_layers_mm[0]-0.2,"Skin",fontsize=ft_txt)
    plt.text(1,thickness_layers_mm[1]-0.2,"Adipose tissue",fontsize=ft_txt)
    plt.text(1,thickness_layers_mm[2]-0.2,"Muscle",fontsize=ft_txt)
    plt.text(1,thickness_layers_mm[2]+0.5,"Placenta",fontsize=ft_txt)

    c = plt.colorbar(im)
    c.set_label("Sensitivity probability (%)",fontsize=ft)
    p.invert_yaxis()
    # plt.xticks(x)
    # plt.yticks(y)
    plt.xlabel("Tissue width (mm)",fontsize=ft)
    plt.ylabel("Tissue depth (mm)",fontsize=ft)
    plt.legend(loc="best",fontsize=ft)

p= plt.subplot(133)

plt.title(title[2],fontsize=ft)
im = plt.contourf(x_interp, y_interp,vec_map[2],cmap='plasma')

plt.plot(x_interp,thickness_layers_mm[0]*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)
plt.plot(x_interp,thickness_layers_mm[1]*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)
plt.plot(x_interp,thickness_layers_mm[2]*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)

plt.plot(10,0,'ko',markersize=12,label="Source")
plt.plot(10+SD_separation_cm_array[SD_separation_id]*10,0,'go',markersize=12,label="Detector")

plt.text(1,thickness_layers_mm[0]-0.2,"Skin",fontsize=ft_txt)
plt.text(1,thickness_layers_mm[1]-0.2,"Adipose tissue",fontsize=ft_txt)
plt.text(1,thickness_layers_mm[2]-0.2,"Muscle",fontsize=ft_txt)
plt.text(1,thickness_layers_mm[2]+0.5,"Placenta",fontsize=ft_txt)

c = plt.colorbar(im)
c.set_label("Sensitivity probability (%)",fontsize=ft)
p.invert_yaxis()
# plt.xticks(x)
# plt.yticks(y)
plt.xlabel("Tissue width (mm)",fontsize=ft)
plt.ylabel("Tissue depth (mm)",fontsize=ft)
plt.legend(loc="best",fontsize=ft)
plt.show()
