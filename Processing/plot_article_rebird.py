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

    # dr_phantom_simu = process_Binning(dr_phantom_simu,binning)

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



## Plot data thickness data and

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

#Size of the vector (integration time, SD separation, Time)
Intensity_measured = []

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
    print(I_interp.shape)


    #calculate SNR
    for binning in np.array([1,2,5,10]):
        binned_signal = process_Binning(I_interp,binning)

        #Get SNR_y
        SNR_y = binned_signal.mean(axis=1) / binned_signal.std(axis=1)


        print("ti ",integration_time_s[j],"binning", binning, "SNR ",SNR_y.mean(), "size signal: ",binned_signal.shape)


# plt.figure()
# for j in range(integration_time_s_Mini_CYRIL.shape[0]):
#
#     file = glob.glob(path+"*"+str(SD_separations_cm_Mini_CYRIL[0])+"cm_"+str(integration_time_s_Mini_CYRIL[j])+"s*")
#     data = scipy.io.loadmat(file[0])
#     plt.plot(Mini_CYRIL_wavelength,data['Spectra'].mean(axis=0),label=str(integration_time_s_Mini_CYRIL[j]))
# plt.legend(loc="best")
# plt.show()

## Load and Plot tissue sensitivity data calculated with RedBird

wavelength = 890

#Load tissue sensitivity
path_data = main_path+"simulations/Redbird/data_article_"+str(wavelength)+"/"
HbT_placenta_array = np.array([15,25,35,50])
SD_separation_cm = np.array([3, 4, 5])

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
        im = ax1.pcolor(T_x,T_y,Placenta_sensitivity_redbird*100, cmap=cmap, vmin=0, vmax=70)

        for x_id,xval in enumerate(x):
            for y_id,yval in enumerate(y):
                ax1.text(xval-0.15,yval,str(int(10000*Placenta_sensitivity_redbird[y_id,x_id])/100),color=colors[i],fontsize=ft)

        cb = fig1.colorbar(im, ax=ax1)
        cb.set_label("Placenta sensitivity (%)",fontsize=ft_label)
        ax1.set_xticks(x)  # Set x ticks to vector values
        ax1.set_yticks(y)  # Set y ticks to vector values
        ax1.set_xlabel("Source detector separation (cm)",fontsize=ft)
        ax1.set_ylabel("Placenta $C_{HbT}$ ($\mu$M)",fontsize=ft)

plt.show()

all_d = np.asarray(all_d)
print(wavelength, int(100*all_d.mean()))

## Calculate the detection probability


wavelength = 780
# Load simulation of phantom measurements
# Load phantom
path = main_path + "simulations/Redbird/"
data_phantom = scipy.io.loadmat(path+'Phantom_Data_'+str(wavelength)+'.mat')
dr_phantom = np.squeeze(data_phantom['Diffuse_reflectance'])


ft = 23
ft_label = 23

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


            val = np.squeeze(data['Diffuse_reflectance'])/Coeff

            for d in range(SD_separation_cm.shape[0]):
                # Detection probability
                sample = np.random.normal(loc = val[d],
                                            scale = sigma_p_phantom_detector[d],
                                            size = (100))
                # Process T-tests
                t,p = scipy.stats.ttest_1samp(sample,mu_p_min[d],
                                            alternative='greater')
                Detection_probability[j,d] = 1-p #Proba to accept alternative hypothesis


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
dr_phantom = np.squeeze(data_phantom['Diffuse_reflectance'])




ft = 23
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
Diffuse_reflectance = data['Diffuse_reflectance']
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

                # Detection probability
                sample = np.random.normal(loc = val,
                                            scale = sigma_p_phantom_detector[d],
                                            size = (100))
                # Process T-tests
                t,p = scipy.stats.ttest_1samp(sample,mu_p_min[d],
                                            alternative='greater')
                Detection_probability[s] = 1-p #Proba to accept alternative hypothesis
                Scanning_probability[s] = (1-p)* S_placenta[d,s]


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
gs = gridspec.GridSpec(3, 1)

ax = fig.add_subplot(gs[0,0])
ax.set_title("Placenta sensitivity",fontsize = ft)
plots = ax.violinplot(np.asarray(S_placenta_m).T*100, showmedians=True, showmeans=False,quantiles=[[0.25,0.75],[0.25,0.75],[0.25,0.75]])
plots['cmedians'].set_colors(c[id])
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
plots = ax.violinplot(np.asarray(P_detection_m).T*100, showmedians=True, showmeans=False)#,quantiles=quantiles_vec)
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)
plots['cmedians'].set_colors(c[id])
# plots['cquantiles'].set_color('r')
ax.set_xticks(np.array([3.5,9.5,15.5]))
ax.set_xticklabels(np.array(["3","4","5"]),fontsize=ft)
ax.set_ylabel("Probability (%)",fontsize = ft)
ax.set_xlabel("Source-detector separation (cm)",fontsize = ft)
ax.axvline(6.5,0,100,color='k',linewidth=3)
ax.axvline(12.5,0,100,color='k',linewidth=3)
ax.grid()
ax.set_ylim(0,100)
ax.legend(handles=patches,loc="best", borderaxespad=0. ,fontsize = ft)


#Scanning probability
ax = fig.add_subplot(gs[2,0])
ax.set_title("Placenta scanning probability",fontsize = ft)
plots = ax.violinplot(np.asarray(P_scanning_m).T*100, showmedians=True, showmeans=False)#,quantiles=quantiles_vec)
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)
plots['cmedians'].set_colors(c[id])
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
dr_phantom = np.squeeze(data_phantom['Diffuse_reflectance'])



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





for i in range(integration_time_s.shape[0]):
    #Plot mu_p_phantom
    mup_min, sigma, binned_signal = get_minimum_sensitivity_threshold_normalization(1,Intensity_measured[i],dr_phantom)

    mup_min_vec.append(mup_min)
    # mes_vec.append(binned_signal)
    # ax1.plot(SD_separations_mm,mup_min,color=colors[id],linewidth=lw,label=device+ " (IT = "+str(integration_time_s[i])+" s)")


    mup_min, sigma, binned_signal = get_minimum_sensitivity_threshold_normalization(nb_binning,Intensity_measured[i],dr_phantom)
    mup_min_vec.append(mup_min)
    # mes_vec.append(binned_signal)
    # ax1.plot(SD_separations_mm,mup_min,color=colors[id+1],linewidth=lw,label=device+ " (IT = "+str(integration_time_s[i])+" s, temporal binning 10)")




    SNR.append(Intensity_measured[i].mean(axis=1)/Intensity_measured[i].std(axis=1))
    mes_vec.append(Intensity_measured[i].std(axis=1))
    label.append("Mini CYRIL (IT = "+str(integration_time_s[i])+" s)")
    id +=1


    sig = process_Binning(Intensity_measured[i],nb_binning)
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

    rects = ax1.bar(xbar + offset, mes_vec[:,i], width, label=label[i])
    # rects = ax1.bar(i,mup_min_vec[0,i], label=label[i])



    # ax_mes.plot(SD_separations_mm,mes_vec[:,i], label=label[i])



ax1.set_title("Minimum detectable diffuse reflectance of Mini CYRIL device",fontsize=ft)
ax1.legend(loc="best",fontsize=ft_legend)
ax1.set_xlabel("Source detector separation (mm)",fontsize=ft)
ax1.set_ylabel("$\phi_{min}$ (a. u.)",fontsize=ft)
ax1.set_xticks(np.arange(SD_separations_mm.shape[0])+5*width/2,np.array(["30","40","50"]),fontsize=ft)
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



##

a = np.random.random(50)
b = np.array([np.sum(a[0:10]),np.sum(a[10:20]),np.sum(a[20:30]),np.sum(a[30:40]),np.sum(a[40:50])])

