import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker, cm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec

import scipy.stats
import glob
from scipy import interpolate
import pandas as pd
import matplotlib.patches as mpatches

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






def get_thickness(percentile):

    # Extract thickness from segmentation files
    info_thickness = info_subject()
    info_thickness = read_thickness_data(info_thickness)
    info_thickness = read_thickness_data(info_thickness,col='.1')


    t_skin = info_thickness.thickness_skin.astype(float)
    t_muscle = info_thickness.thickness_muscle.astype(float)
    t_adipose = info_thickness.thickness_adipose.astype(float)
    Fitzpatrick = info_thickness.Fitzpatrick_scale.astype(float)



    id_nan = ~ np.isnan(t_skin)
    t_skin = t_skin[id_nan]
    t_adipose = t_adipose[id_nan]
    t_muscle = t_muscle[id_nan]
    Fitzpatrick = Fitzpatrick[id_nan]

    id_nan = ~ np.isnan(t_muscle)
    t_skin = t_skin[id_nan]
    t_adipose = t_adipose[id_nan]
    t_muscle = t_muscle[id_nan]
    Fitzpatrick = Fitzpatrick[id_nan]
    fmel = convert_Fitzpatrick_scale_f_mel(Fitzpatrick)


    config = display_config()
    config.id_skin = np.percentile(t_skin*10,percentile)
    config.id_adipose = np.percentile(t_adipose*10,percentile)
    config.id_muscle = np.percentile(t_muscle*10,percentile)
    config.id_mel = np.percentile(fmel,percentile)

    return config



def get_thickness_indexes(percentile, thickness_skin_mm = -1):

    # Extract thickness from segmentation files
    info_thickness = info_subject()
    info_thickness = read_thickness_data(info_thickness)
    info_thickness = read_thickness_data(info_thickness,col='.1')


    t_skin = info_thickness.thickness_skin.astype(float)
    t_muscle = info_thickness.thickness_muscle.astype(float)
    t_adipose = info_thickness.thickness_adipose.astype(float)
    Fitzpatrick = info_thickness.Fitzpatrick_scale.astype(float)



    id_nan = ~ np.isnan(t_skin)
    t_skin = t_skin[id_nan]
    t_adipose = t_adipose[id_nan]
    t_muscle = t_muscle[id_nan]
    Fitzpatrick = Fitzpatrick[id_nan]

    id_nan = ~ np.isnan(t_muscle)
    t_skin = t_skin[id_nan]
    t_adipose = t_adipose[id_nan]
    t_muscle = t_muscle[id_nan]
    Fitzpatrick = Fitzpatrick[id_nan]
    fmel = convert_Fitzpatrick_scale_f_mel(Fitzpatrick)


    config = display_config()

    if thickness_skin_mm == -1:
        config.id_skin = np.argmin(np.abs(skin_thickness_array - np.percentile(t_skin*10,percentile)))
    else:
        config.id_skin = np.argmin(np.abs(skin_thickness_array - thickness_skin_mm))

    config.id_adipose = np.argmin(np.abs(adipose_thickness_array - np.percentile(t_adipose*10,percentile)))
    config.id_muscle = np.argmin(np.abs(muscle_thickness_array - np.percentile(t_muscle*10,percentile)))
    config.id_mel = np.argmin(np.abs(f_melanosome - np.percentile(fmel,percentile)))

    return config


def plot_detection_proba(title,Proba,Sensor_signal,SD_separations_mm, conf = display_config(), mode='Adipose vs Muscle',mu_p_min=0,noise_level=0,ft_label= 16,ft_title = 18):


    #Colormaps
    cmap = cm.plasma
    # cmap = cm.jet

    alpha = 1



    #angle for vizualization
    elev = 20
    azim = 60
    roll = 0

    plt.rcParams.update({'font.size': ft_label})


    #Mode adipose vs muscle
    if(mode == 'Adipose vs Muscle'):
        conf.id_adipose = np.arange(adipose_thickness_array.shape[0])
        conf.id_muscle = np.arange(muscle_thickness_array.shape[0])

        # title += " - $HbT$ muscle "+str(HbT_muscle_array[conf.id_HbT_m])+ " $\mu Mol$ - $HbT$ placenta "+str(HbT_placenta_array[conf.id_HbT_p])+ " $\mu Mol$"+"\n$SatO_2$ muscle "+str(int(SatO2_muscle_array[conf.id_sat_m]*100))+ "% - $SatO_2$ placenta "+str(int(SatO2_placenta_array[conf.id_sat_p]*100))+ "% - Thickness skin "+str(skin_thickness_array[conf.id_skin]) + " mm - Fitzpatrick skin type "+ convert_f_mel_Fitzpatrick_scale(f_melanosome[conf.id_mel])
        title += "\nInfluence of adipose tissue and muscle thickness"

        #meshgrid(y,x)
        x = adipose_thickness_array
        y = muscle_thickness_array


        #labels
        x_label = "Adipose thickness (mm)"
        y_label = "Muscle thickness (mm)"

        #Select proba
        P_detection = Proba[conf.id_skin, :,:, conf.id_sat_m, conf.id_sat_p, conf.id_mel, conf.id_HbT_m, conf.id_HbT_p, :]
        Signal = Sensor_signal[conf.id_skin, :,:, conf.id_sat_m, conf.id_sat_p, conf.id_mel, conf.id_HbT_m, conf.id_HbT_p, :]







    #Mode HbT
    if(mode == "HbT"):
        conf.id_HbT_m = np.arange(HbT_muscle_array.shape[0])
        conf.id_HbT_p = np.arange(HbT_placenta_array.shape[0])

        # title += " - Thickness skin "+str(skin_thickness_array[conf.id_skin]) + " mm - Thickness adipose "+ str(adipose_thickness_array[conf.id_adipose])+ " mm\nThickness muscle "+str(muscle_thickness_array[conf.id_muscle])+" mm - Fitzpatrick skin type "+ convert_f_mel_Fitzpatrick_scale(f_melanosome[conf.id_mel]) + " - $SatO_2$ muscle "+str(int(SatO2_muscle_array[conf.id_sat_m]*100))+ "% - $SatO_2$ placenta "+str(int(SatO2_placenta_array[conf.id_sat_p]*100))

        title += "\nInfluence of muscle and placenta blood volume"

        #meshgrid(y,x)
        x = HbT_muscle_array
        y = HbT_placenta_array

        #labels
        x_label = "$HbT$ muscle (µM)"
        y_label = "$HbT$ placenta (µM)"


        #Select proba
        P_detection = Proba[conf.id_skin, conf.id_adipose, conf.id_muscle, conf.id_sat_m, conf.id_sat_p, conf.id_mel, :, :, :]
        Signal = Sensor_signal[conf.id_skin, conf.id_adipose, conf.id_muscle, conf.id_sat_m, conf.id_sat_p, conf.id_mel, :, :, :]





    #Mode skin
    if(mode == "Skin"):
        conf.id_skin = np.arange(skin_thickness_array.shape[0])
        conf.id_mel = np.arange(f_melanosome.shape[0])

        # title += " - Thickness adipose "+ str(adipose_thickness_array[conf.id_adipose])+ " mm - Thickness muscle "+str(muscle_thickness_array[conf.id_muscle])+" mm\n$HbT$ muscle "+str(HbT_muscle_array[conf.id_HbT_m])+ " $\mu Mol$ - $HbT$ placenta "+str(HbT_placenta_array[conf.id_HbT_p])+ " $\mu Mol$ - $SatO_2$ muscle "+str(int(SatO2_muscle_array[conf.id_sat_m]*100))+ "% - $SatO_2$ placenta "+str(int(SatO2_placenta_array[conf.id_sat_p]*100))+"%"

        title += "\nInfluence of skin tone and thickness"

        #meshgrid(y,x)
        x = skin_thickness_array
        y = f_melanosome


        #labels
        x_label = "Skin thickness (mm)"
        y_label = "Melanosome volume fraction (%)"


        #Select proba
        P_detection = Proba[:, conf.id_adipose, conf.id_muscle, conf.id_sat_m, conf.id_sat_p, :, conf.id_HbT_m, conf.id_HbT_p, :]
        Signal = Sensor_signal[:, conf.id_adipose, conf.id_muscle, conf.id_sat_m, conf.id_sat_p, :, conf.id_HbT_m, conf.id_HbT_p, :]




    #Meshgrid
    x_interp = np.linspace(x[0],x[-1],100)
    y_interp = np.linspace(y[0],y[-1],100)

    T_yy,T_xx = np.meshgrid(y_interp,x_interp)
    T_y,T_x = np.meshgrid(y,x)



    # Plot
    nb_row = 2
    nb_col = SD_separations_mm.shape[0]



    fig = plt.figure()
    plt.suptitle(title,fontsize=ft_title)

    for i in range(SD_separations_mm.shape[0]):
        #Probability
        interp = interpolate.RegularGridInterpolator((x, y),P_detection[:,:,i],bounds_error=False, fill_value=None)
        P_interp = interp((T_xx,T_yy))
        P_interp[P_interp<0]=0
        P_interp[P_interp>100]=100

        plt.subplot(nb_row,nb_col,i+1)
        plt.title("Detector "+str(int(SD_separations_mm[i]))+" mm\nDetection probability",fontsize=ft_title)
        im = plt.pcolor(T_xx,T_yy,P_interp, vmin=0, vmax=100, cmap=cmap)
        plt.xlabel(x_label,fontsize=ft_label)
        plt.ylabel(y_label,fontsize=ft_label)
        cb = plt.colorbar(im)
        cb.set_label("Probability (%)",fontsize=ft_label)

        #Sensor signal
        interp = interpolate.RegularGridInterpolator((x, y),Signal[:,:,i],bounds_error=False, fill_value=None)
        Signal_interp = interp((T_xx,T_yy))

        plt.subplot(nb_row,nb_col,SD_separations_mm.shape[0]+i+1)
        plt.title("Detector "+str(int(SD_separations_mm[i]))+" mm\nDiffuse reflectance",fontsize=ft_title)

        vmax = Signal_interp.max()
        if vmax<mu_p_min[i]:
            vmax = 1.05*mu_p_min[i]
        im = plt.pcolor(T_xx,T_yy,Signal_interp, cmap=cmap,vmin=0,vmax=vmax)

        plt.xlabel(x_label,fontsize=ft_label)
        plt.ylabel(y_label,fontsize=ft_label)
        cb = plt.colorbar(im)
        cb.ax.axhline(mu_p_min[i], c='w')
        cb.ax.axhline(noise_level[i], c='k')
        cb.set_label("$\phi$ (a. u.)",fontsize=ft_label)

    plt.show()


def plot_Tissue_sensitivity(Proba, title, det_mm, conf = display_config(), mode='Adipose vs Muscle',mode_display='',levels = [0,0,0,0], ft_label=12, ft_title=14):


    #Init title
    title += " - Source/Detector separation "+str(det_mm)+"mm"


    #Colormaps
    # cmap = cm.plasma
    cmap = cm.jet
    alpha = 1

    #Detector
    SD_separations_mm = np.arange(10,90,10)
    id_det = np.where(SD_separations_mm==det_mm)[0].item()


    #angle for vizualization
    elev = 20
    azim = 60
    roll = 0

    plt.rcParams.update({'font.size': ft_label})


    #Mode adipose vs muscle
    if(mode == 'Adipose vs Muscle'):
        conf.id_adipose = np.arange(adipose_thickness_array.shape[0])
        conf.id_muscle = np.arange(muscle_thickness_array.shape[0])

        #title += "\n$HbT$ muscle "+str(HbT_muscle_array[conf.id_HbT_m])+ "$\mu Mol$ - $HbT$ placenta "+str(HbT_placenta_array[conf.id_HbT_p])+ "$\mu Mol$"+"\n$SatO_2$ muscle "+str(int(SatO2_muscle_array[conf.id_sat_m]*100))+ "% - $SatO_2$ placenta "+str(int(SatO2_placenta_array[conf.id_sat_p]*100))+ "%\n Thickness skin "+str(skin_thickness_array[conf.id_skin]) + " mm - Fitzpatrick skin type "+ convert_f_mel_Fitzpatrick_scale(f_melanosome[conf.id_mel])

        title += "\nInfluence of adipose tissue and muscle thickness"

        #meshgrid(y,x)
        x = adipose_thickness_array
        y = muscle_thickness_array

        #labels
        x_label = "Adipose thickness (mm)"
        y_label = "Muscle thickness (mm)"

        #Select proba
        P_skin = 100*Proba[conf.id_skin, :,:, conf.id_sat_m, conf.id_sat_p, conf.id_mel, conf.id_HbT_m, conf.id_HbT_p, id_det, 0]

        P_adipose = 100*Proba[conf.id_skin, :,:, conf.id_sat_m, conf.id_sat_p, conf.id_mel, conf.id_HbT_m, conf.id_HbT_p, id_det, 1]

        P_muscle = 100*Proba[conf.id_skin, :,:, conf.id_sat_m, conf.id_sat_p, conf.id_mel, conf.id_HbT_m, conf.id_HbT_p, id_det, 2]

        P_placenta = 100*Proba[conf.id_skin, :,:, conf.id_sat_m, conf.id_sat_p, conf.id_mel, conf.id_HbT_m, conf.id_HbT_p, id_det, 3]







    #Mode HbT
    if(mode == "HbT"):
        conf.id_HbT_m = np.arange(HbT_muscle_array.shape[0])
        conf.id_HbT_p = np.arange(HbT_placenta_array.shape[0])

        # title += "\nThickness skin "+str(skin_thickness_array[conf.id_skin]) + "mm - Thickness adipose "+ str(adipose_thickness_array[conf.id_adipose])+ " mm" + "\n Thickness muscle "+str(muscle_thickness_array[conf.id_muscle])+" mm - Fitzpatrick skin type "+ convert_f_mel_Fitzpatrick_scale(f_melanosome[conf.id_mel]) + "\n$SatO_2$ muscle "+str(int(SatO2_muscle_array[conf.id_sat_m]*100))+ "% - $SatO_2$ placenta "+str(int(SatO2_placenta_array[conf.id_sat_p]*100))

        title += "\nInfluence of muscle and placenta blood volume"

        #meshgrid(y,x)
        x = HbT_muscle_array
        y = HbT_placenta_array


        #labels
        x_label = "$HbT$ in muscle (µM)"
        y_label = "$HbT$ in placenta (µM)"


        #Select proba
        P_skin = 100*Proba[conf.id_skin, conf.id_adipose, conf.id_muscle, conf.id_sat_m, conf.id_sat_p, conf.id_mel, :,:, id_det, 0]

        P_adipose = 100*Proba[conf.id_skin, conf.id_adipose, conf.id_muscle, conf.id_sat_m, conf.id_sat_p, conf.id_mel, :,:, id_det, 1]

        P_muscle = 100*Proba[conf.id_skin, conf.id_adipose, conf.id_muscle, conf.id_sat_m, conf.id_sat_p, conf.id_mel, :,:, id_det, 2]

        P_placenta = 100*Proba[conf.id_skin, conf.id_adipose, conf.id_muscle, conf.id_sat_m, conf.id_sat_p, conf.id_mel, :,:, id_det, 3]






    #Mode skin
    if(mode == "Skin"):
        conf.id_skin = np.arange(skin_thickness_array.shape[0])
        conf.id_mel = np.arange(f_melanosome.shape[0])

        # title += "\nThickness adipose "+ str(adipose_thickness_array[conf.id_adipose])+ " mm" + "\n Thickness muscle "+str(muscle_thickness_array[conf.id_muscle])+" mm" + "\n$HbT$ muscle "+str(HbT_muscle_array[conf.id_HbT_m])+ "$\mu Mol$ - $HbT$ placenta "+str(HbT_placenta_array[conf.id_HbT_p])+ "$\mu Mol$"+"\n$SatO_2$ muscle "+str(int(SatO2_muscle_array[conf.id_sat_m]*100))+ "% - $SatO_2$ placenta "+str(int(SatO2_placenta_array[conf.id_sat_p]*100))+"%"

        title += "\nInfluence of skin tone and thickness"

        #meshgrid(y,x)
        x = skin_thickness_array
        y = f_melanosome


        #labels
        x_label = "Skin thickness (mm)"
        y_label = "Melanosome volume fraction (%)"


        #Select proba
        P_skin = 100*Proba[:, conf.id_adipose, conf.id_muscle, conf.id_sat_m, conf.id_sat_p, :, conf.id_HbT_m, conf.id_HbT_p, id_det, 0]

        P_adipose = 100*Proba[:, conf.id_adipose, conf.id_muscle, conf.id_sat_m, conf.id_sat_p, :, conf.id_HbT_m, conf.id_HbT_p, id_det, 1]

        P_muscle = 100*Proba[:, conf.id_adipose, conf.id_muscle, conf.id_sat_m, conf.id_sat_p, :, conf.id_HbT_m, conf.id_HbT_p, id_det, 2]

        P_placenta = 100*Proba[:, conf.id_adipose, conf.id_muscle, conf.id_sat_m, conf.id_sat_p, :, conf.id_HbT_m, conf.id_HbT_p, id_det, 3]


























    #Interpolation
    Nb_interp = 10
    x_interp = np.linspace(x[0],x[-1],Nb_interp)
    y_interp = np.linspace(y[0],y[-1],Nb_interp)
    T_y,T_x = np.meshgrid(y_interp,x_interp)


    interp_skin = interpolate.RegularGridInterpolator((x, y),P_skin,bounds_error=False, fill_value=None)
    interp_adipose = interpolate.RegularGridInterpolator((x, y),P_adipose,bounds_error=False, fill_value=None)
    interp_muscle = interpolate.RegularGridInterpolator((x, y),P_muscle,bounds_error=False, fill_value=None)
    interp_placenta = interpolate.RegularGridInterpolator((x, y),P_placenta,bounds_error=False, fill_value=None)


    P_skin = interp_skin((T_x,T_y))
    P_adipose = interp_adipose((T_x,T_y))
    P_muscle = interp_muscle((T_x,T_y))
    P_placenta = interp_placenta((T_x,T_y))






    # Plot the surface

    if(mode_display == '3d'):
        fig = plt.figure()
        plt.suptitle(title,fontsize=ft_title)

        ax = fig.add_subplot(2,2,1,projection='3d')
        ax.set_title("Skin",fontsize=ft_title)
        surf = ax.plot_surface(T_x, T_y, P_skin, cmap=cmap, linewidth=0, antialiased=True,vmin = 0,vmax= 100,alpha=alpha)

        # Customize the z axis.
        ax.set_zlim(0, 100)
        ax.set_xlabel(x_label,fontsize=ft_label)
        ax.set_ylabel(y_label,fontsize=ft_label)
        ax.view_init(elev, azim, roll)

        ax = fig.add_subplot(2,2,2,projection='3d')
        ax.set_title("Adipose tissue",fontsize=ft_title)
        surf = ax.plot_surface(T_x, T_y, P_adipose, cmap=cmap, linewidth=0, antialiased=True,vmin = 0,vmax= 100,alpha=alpha)

        # Customize the z axis.
        ax.set_zlim(0, 100)
        ax.set_xlabel(x_label,fontsize=ft_label)
        ax.set_ylabel(y_label,fontsize=ft_label)
        ax.view_init(elev, azim, roll)


        ax = fig.add_subplot(2,2,3,projection='3d')
        ax.set_title("Muscle",fontsize=ft_title)
        surf = ax.plot_surface(T_x, T_y, P_muscle, cmap=cmap, linewidth=0, antialiased=True,vmin = 0,vmax= 100,alpha=alpha)

        # Customize the z axis.
        ax.set_zlim(0, 100)
        ax.set_xlabel(x_label,fontsize=ft_label)
        ax.set_ylabel(y_label,fontsize=ft_label)
        ax.view_init(elev, azim, roll)

        ax = fig.add_subplot(2,2,4,projection='3d')
        ax.set_title("Placenta",fontsize=ft_title)
        surf = ax.plot_surface(T_x, T_y, P_placenta, cmap=cmap, linewidth=0, antialiased=True,vmin = 0,vmax= 100,alpha=alpha)

        # Customize the z axis.
        ax.set_zlim(0, 100)
        ax.set_xlabel(x_label,fontsize=ft_label)
        ax.set_ylabel(y_label,fontsize=ft_label)
        ax.view_init(elev, azim, roll)

        plt.show()

    else:


        vmin = np.min([P_skin.min(),P_adipose.min(),P_muscle.min(),P_placenta.min()])
        vmax = np.max([P_skin.max(),P_adipose.max(),P_muscle.max(),P_placenta.max()])
        # vmin = 0

        fig = plt.figure()
        plt.suptitle(title,fontsize=ft_title)

        plt.subplot(221)
        plt.title("Skin",fontsize=ft_title)
        if np.size(levels[0]) == 1:
            im = plt.contourf(T_x, T_y,P_skin,cmap='plasma')
        else:
            im = plt.contourf(T_x, T_y,P_skin,cmap='plasma',levels=levels[0])

        plt.xlabel(x_label,fontsize=ft_label)
        plt.ylabel(y_label,fontsize=ft_label)
        c=plt.colorbar(im)
        c.set_label("Probability (%)")

        plt.subplot(222)
        plt.title("Adipose tissue",fontsize=ft_title)
        if np.size(levels[1]) == 1:
            im = plt.contourf(T_x, T_y,P_adipose,cmap='plasma')
        else:
            im = plt.contourf(T_x, T_y,P_adipose,cmap='plasma',levels=levels[1])
        plt.xlabel(x_label,fontsize=ft_label)
        plt.ylabel(y_label,fontsize=ft_label)
        c=plt.colorbar(im)
        c.set_label("Probability (%)")

        plt.subplot(223)
        plt.title("Muscle",fontsize=ft_title)
        if np.size(levels[2]) == 1:
            im = plt.contourf(T_x, T_y,P_muscle,cmap='plasma')
        else:
            im = plt.contourf(T_x, T_y,P_muscle,cmap='plasma',levels=levels[2])
        plt.xlabel(x_label,fontsize=ft_label)
        plt.ylabel(y_label,fontsize=ft_label)
        c=plt.colorbar(im)
        c.set_label("Probability (%)")


        plt.subplot(224)
        plt.title("Placenta",fontsize=ft_title)
        if np.size(levels[3]) == 1:
            im = plt.contourf(T_x, T_y,P_placenta,cmap='plasma')
        else:
            im = plt.contourf(T_x, T_y,P_placenta,cmap='plasma',levels=levels[3])
        plt.xlabel(x_label,fontsize=ft_label)
        plt.ylabel(y_label,fontsize=ft_label)
        c=plt.colorbar(im)
        c.set_label("Probability (%)")

        plt.show()




def calculate_mua(w,W,F,C_HbT,SatO2):

    # unit in cm-1/Mol
    path = main_path+"/simulations/spectra/"
    wavelength = np.loadtxt(path+"lambda.txt")
    eps_Hb = np.loadtxt(path+"eps_Hb.txt")
    eps_HbO2 = np.loadtxt(path+"eps_HbO2.txt")

    mua_W = np.loadtxt(path+"mua_H2O.txt")
    mua_F = np.loadtxt(path+"mua_Fat.txt")

    _eps_Hb = scipy.interpolate.interp1d(wavelength,eps_Hb, kind='cubic')(w)
    _eps_HbO2 = scipy.interpolate.interp1d(wavelength,eps_HbO2, kind='cubic')(w)
    _mua_W = scipy.interpolate.interp1d(wavelength,mua_W, kind='cubic')(w)
    _mua_F = scipy.interpolate.interp1d(wavelength,mua_F, kind='cubic')(w)

    #Calculate absorption coefficient
    mua_placenta = W*_mua_W + F*_mua_F + np.log(10)*C_HbT*(1-SatO2)*_eps_Hb + np.log(10)*C_HbT*SatO2*_eps_HbO2 #cm-1

    mua_placenta = 0.1*mua_placenta # convert into mm-1

    return mua_placenta



def calculate_scanning_proba(val_skin_mm, val_adipose_mm, val_muscle_mm, device_name, binning, int_time_s):

    #load tissue sensitivity
    data = np.load(main_path+"simulations/tissue_sensitivity.npz")
    Tissue_sensitivity_indexes = data['Tissue_sensitivity_indexes']
    skin_thickness_array = data['skin_thickness_array']
    adipose_thickness_array = data['adipose_thickness_array']
    muscle_thickness_array = data['muscle_thickness_array']
    SatO2_muscle_array = data['SatO2_muscle_array']
    SatO2_placenta_array = data['SatO2_placenta_array']
    f_melanosome = data['f_melanosome']
    HbT_muscle_array = data['HbT_muscle_array']
    HbT_placenta_array = data['HbT_placenta_array']
    SD_separation_cm = data['SD_separation_cm']


    id_skin = np.argmin(np.abs(skin_thickness_array - val_skin_mm))
    id_adipose = np.argmin(np.abs(adipose_thickness_array - val_adipose_mm))
    id_muscle = np.argmin(np.abs(muscle_thickness_array - val_muscle_mm))

    Scanning_proba_muscle = []
    Scanning_proba_placenta = []
    labels = []

    #Load detection proba
    file = main_path+"simulations/detection_proba/detection_probability_"+device_name+"_TI_"+str(int_time_s)+"s_binning_"+str(binning)+"_camera_model_0.npz"
    data = np.load(file)
    Detection_probability = data['Detection_probability']


    if device_name == "Mini CYRIL":
        id_detector_det = np.array([0,1,2])
        id_detector_sensi = np.array([2,3,4])
    if device_name == "CYRIL":
        id_detector_det = np.array([0,1,2,3,4])
        id_detector_sensi = np.array([0,1,2,3,4])


    Scanning_proba_muscle = Tissue_sensitivity_indexes[id_skin, id_adipose, id_muscle, :, :, :,:,:,id_detector_sensi, 2]*Detection_probability[id_skin, id_adipose, id_muscle, :, :, :,:,:,id_detector_det]

    Scanning_proba_placenta = Tissue_sensitivity_indexes[id_skin, id_adipose, id_muscle, :, :, :,:,:,id_detector_sensi, 3]*Detection_probability[id_skin, id_adipose, id_muscle, :, :, :,:,:,id_detector_det]


    return Scanning_proba_muscle, Scanning_proba_placenta


## Plots thickness distribution
ft = 16


# device_name = "Mini CYRIL"
device_name = "CYRIL"


title_distrib = "sensitivity indexes"
# title_distrib = "scanning probability"




# Extract thickness from segmentation files
info_thickness = info_subject()
info_thickness = read_thickness_data(info_thickness)
info_thickness = read_thickness_data(info_thickness,col='.1')
info_thickness = read_thickness_data(info_thickness,col='.2')


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


##

info = []
info.append(t_skin)
info.append(t_adipose)
info.append(t_muscle)

ft = 22
plt.close('all')
plt.figure()
labels = np.array(["Skin","Adipose tissue", "Muscle"])
quantiles = [[0.25,0.70],[0.25,0.70],[0.25,0.70]]

plt.subplot(221)
plots = plt.violinplot(info,showmedians=True, showmeans=False,quantiles=quantiles)
plots['cquantiles'].set_color('k')
# plots['cmeans'].set_color('r')
plots['cmedians'].set_color('r')
plt.xticks(np.arange(1,len(info)+1),labels,fontsize = ft)
plt.ylabel("Thickness (cm)",fontsize = ft)
plt.grid()
plt.title("Tissue thickness",fontsize = ft)
plt.yticks(fontsize=ft)


plt.subplot(222)
plt.title("Distance between skin and placenta",fontsize = ft)
plots = plt.violinplot(dist_to_placenta,showmedians=True, showmeans=False,quantiles=[0.25,0.75])
plots['cquantiles'].set_color('k')
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
plots = plt.violinplot(Gestation_subjects,showmedians=True, showmeans=False,quantiles=[0.25,0.75])
plots['cquantiles'].set_color('k')
# plots['cmeans'].set_color('r')
plots['cmedians'].set_color('r')
plt.ylabel("Gestational age (weeks)",fontsize = ft)
plt.xticks([],'')
plt.grid()
plt.yticks(fontsize=ft)
plt.show()


## Load and Plot tissue sensitivity data

device_name = "CYRIL"
# device_name = "Mini CYRIL"

#Save tissue sensitivity
data = np.load(main_path+"simulations/tissue_sensitivity.npz")



Tissue_sensitivity_indexes = data['Tissue_sensitivity_indexes']
skin_thickness_array = data['skin_thickness_array']
adipose_thickness_array = data['adipose_thickness_array']
muscle_thickness_array = data['muscle_thickness_array']
SatO2_muscle_array = data['SatO2_muscle_array']
SatO2_placenta_array = data['SatO2_placenta_array']
f_melanosome = data['f_melanosome']
HbT_muscle_array = data['HbT_muscle_array']
HbT_placenta_array = data['HbT_placenta_array']
SD_separation_cm = data['SD_separation_cm']



# Plot tissue sensitivity


#Get median thicknesses and fmel config
config = get_thickness_indexes(50)

config.id_HbT_m = 1
config.id_HbT_p = 1
config.id_sat_p = 2
config.id_sat_m = 1

det_mm = 30

mode_display = ''

# mode = 'Adipose vs Muscle'
mode = 'HbT'
# mode = 'Skin'

#set levels
if mode == 'Adipose vs Muscle':
    levels = [] #0 is auto level
    levels.append(0) #skin
    levels.append(np.arange(0,110,10)) #Adipose tissue
    levels.append(np.arange(0,110,10)) #Muscle
    levels.append(np.arange(0,110,10)) #Placenta
if mode == 'HbT':
    levels = [0,0,0,0]

if mode == 'Skin':
    levels = [0,0,0,0]


plt.close('all')
ft_label= 16
ft_title = 18
plot_Tissue_sensitivity(Tissue_sensitivity_indexes, "Tissue sensitivity", det_mm, conf = config, mode = mode,mode_display = mode_display, levels = levels, ft_label=ft_label, ft_title=ft_title)


## Load and Plot tissue sensitivity data v2

device_name = "CYRIL"
# device_name = "Mini CYRIL"

#Save tissue sensitivity
data = np.load(main_path+"simulations/tissue_sensitivity.npz")

ft = 17
ft_label = 16

Tissue_sensitivity_indexes = data['Tissue_sensitivity_indexes']
skin_thickness_array = data['skin_thickness_array']
adipose_thickness_array = data['adipose_thickness_array']
muscle_thickness_array = data['muscle_thickness_array']
SatO2_muscle_array = data['SatO2_muscle_array']
SatO2_placenta_array = data['SatO2_placenta_array']
f_melanosome = data['f_melanosome']
HbT_muscle_array = data['HbT_muscle_array']
HbT_placenta_array = data['HbT_placenta_array']
SD_separation_cm = data['SD_separation_cm']

idx_SD = np.array([2,3,4])
Placenta_sensitivity = Tissue_sensitivity_indexes[:,:,:,:,:,:,:,:,idx_SD,-1]
SD_separation_cm = SD_separation_cm[idx_SD]



skin_thickness_subject_mm = np.array([1, 2, 3])
adipose_thickness_subject_mm = np.array([2, 4, 4])
muscle_thickness_subject_mm = np.array([7, 10, 13])

dist_to_plancenta_subjects_mm = skin_thickness_subject_mm + adipose_thickness_subject_mm + muscle_thickness_subject_mm



plt.close('all')
fig1 = plt.figure()
plt.suptitle("Placenta sensitivity",fontsize=ft)
# fig2 = plt.figure()

cmap = cm.plasma

colors = ['k','w','w','w']

for id_mel in range(3):
    for i in range(muscle_thickness_subject_mm.shape[0]):


        #Display  (HbT placenta,detectors)
        config = display_config()
        config.id_skin = np.argmin(np.abs(skin_thickness_array-skin_thickness_subject_mm[i]))
        config.id_adipose = np.argmin(np.abs(adipose_thickness_array-adipose_thickness_subject_mm[i]))
        config.id_muscle = np.argmin(np.abs(muscle_thickness_array-muscle_thickness_subject_mm[i]))

        config.id_HbT_m = 1
        config.id_sat_p = 2
        config.id_sat_m = 1
        config.id_mel = id_mel


        #shape (fmel,detectors)
        proba = Placenta_sensitivity[config.id_skin, config.id_adipose, config.id_muscle, config.id_sat_m, config.id_sat_p, config.id_mel, config.id_HbT_m, :, :]


        y = HbT_placenta_array
        x = SD_separation_cm.copy()
        T_x,T_y = np.meshgrid(x,y)

        ax1 = fig1.add_subplot(3,3,i+1+(3*(id_mel)))
        ax1.set_title("Placenta depth "+str(dist_to_plancenta_subjects_mm[i])+" mm",fontsize=ft)
        im = ax1.pcolor(T_x,T_y,100*proba, cmap=cmap, vmin=0, vmax=50)

        for x_id,xval in enumerate(x):
            for y_id,yval in enumerate(y):
                ax1.text(xval-0.15,yval,str(int(10000*proba[y_id,x_id])/100),color=colors[i],fontsize=ft)

        cb = fig1.colorbar(im, ax=ax1)
        cb.set_label("Placenta sensitivity (%)",fontsize=ft_label)
        ax1.set_xticks(x)  # Set x ticks to vector values
        ax1.set_yticks(y)  # Set y ticks to vector values
        ax1.set_xlabel("Source detector separation (cm)",fontsize=ft)
        ax1.set_ylabel("Placenta blood volume ($\mu$M)",fontsize=ft)

plt.show()



## Load and Plot detection probabilities v2
plt.close('all')

# device_name = "CYRIL"
device_name = "Mini CYRIL"

ft_label= 16
ft_title = 18
camera_model = 0

# mode_array = np.array(['Adipose vs Muscle','HbT','Skin'])
# mode_array = np.array(['Skin'])
int_time_s = 1
binning = 1

title = "Detection probability - Integration time "+str(int_time_s)+" binning "+str(binning)



#Load detection proba
data = np.load(main_path+"simulations/detection_proba/detection_probability_"+device_name+"_TI_"+str(int_time_s)+"s_binning_"+str(binning)+"_camera_model_"+str(camera_model)+".npz")


Detection_probability = 100*data['Detection_probability'] #convert in %
MC_noise_proba = 100*data['MC_noise_proba']
detector_signal = data['detector_signal']
Diffuse_reflectance = data['Diffuse_reglectance']
mu_p_min = data['mu_p_min']
noise_level = data['noise_level']

skin_thickness_array = data['skin_thickness_array']
adipose_thickness_array = data['adipose_thickness_array']
muscle_thickness_array = data['muscle_thickness_array']
SatO2_muscle_array = data['SatO2_muscle_array']
SatO2_placenta_array = data['SatO2_placenta_array']
f_melanosome = data['f_melanosome']
HbT_muscle_array = data['HbT_muscle_array']
HbT_placenta_array = data['HbT_placenta_array']
SD_separation_cm = data['SD_separation_cm']


#Multiplyt MC proba to diffuse reflectance and detection proba
# Detection_probability = np.multiply(Detection_probability,MC_noise_proba)
# detector_signal = np.multiply(detector_signal,MC_noise_proba)


#Choose detector idx
det_idx = np.arange(SD_separation_cm.shape[0])



SD_separation_cm = SD_separation_cm[det_idx]
mu_p_min = mu_p_min[det_idx]
noise_level = noise_level[det_idx]
detector_signal = detector_signal[:,:,:,:,:,:,:,:,det_idx]
Detection_probability = Detection_probability[:,:,:,:,:,:,:,:,det_idx]
MC_noise_proba = MC_noise_proba[:,:,:,:,:,:,:,:,det_idx]

# Detection_probability[Detection_probability<95]=0
# dist_to_plancenta_subjects_mm = np.array([13, 16, 21])


skin_thickness_subject_mm = np.array([1, 2, 3])
adipose_thickness_subject_mm = np.array([2, 4, 4])
muscle_thickness_subject_mm = np.array([7, 10, 13])

dist_to_plancenta_subjects_mm = skin_thickness_subject_mm + adipose_thickness_subject_mm + muscle_thickness_subject_mm
print(dist_to_plancenta_subjects_mm)




ft = 18
cmap = cm.plasma


plt.close('all')
fig1 = plt.figure()
plt.suptitle(title,fontsize=ft)


fig2 = plt.figure()

fig3 = plt.figure()


for id_mel in range(3):
    for i in range(muscle_thickness_subject_mm.shape[0]):


        #Display  (HbT placenta,detectors)
        config = display_config()
        config.id_skin = np.argmin(np.abs(skin_thickness_array-skin_thickness_subject_mm[i]))
        config.id_adipose = np.argmin(np.abs(adipose_thickness_array-adipose_thickness_subject_mm[i]))
        config.id_muscle = np.argmin(np.abs(muscle_thickness_array-muscle_thickness_subject_mm[i]))

        config.id_HbT_m = 1
        config.id_sat_p = 2
        config.id_sat_m = 1
        config.id_mel = id_mel


        #shape (fmel,detectors)
        P_detection = Detection_probability[config.id_skin, config.id_adipose, config.id_muscle, config.id_sat_m, config.id_sat_p, config.id_mel, config.id_HbT_m, :, :]
        I = detector_signal[config.id_skin, config.id_adipose, config.id_muscle, config.id_sat_m, config.id_sat_p, config.id_mel, config.id_HbT_m, :, :]
        P_noise = MC_noise_proba[config.id_skin, config.id_adipose, config.id_muscle, config.id_sat_m, config.id_sat_p, config.id_mel, config.id_HbT_m, :, :]

        y = HbT_placenta_array
        x = SD_separation_cm.copy()
        T_x,T_y = np.meshgrid(x,y)

        ax1 = fig1.add_subplot(3,3,i+1+(3*(id_mel)))
        ax1.set_title("Placenta depth "+str(dist_to_plancenta_subjects_mm[i])+" mm",fontsize=ft)
        im = ax1.pcolor(T_x,T_y,P_detection, vmin=0, vmax=100, cmap=cmap)
        cb = fig1.colorbar(im, ax=ax1)
        cb.set_label("Detection probability (%)",fontsize=ft_label)
        ax1.set_xticks(x)  # Set x ticks to vector values
        ax1.set_yticks(y)  # Set y ticks to vector values
        ax1.set_xlabel("Source detector separation (cm)",fontsize=ft)
        ax1.set_ylabel("Placenta blood volume ($\mu$M)",fontsize=ft)


        for x_id,xval in enumerate(x):
            for y_id,yval in enumerate(y):
                if P_detection[y_id,x_id]<50:
                    col = 'w'
                    offset = 0
                else:
                    col = 'k'
                    offset = 0.1
                ax1.text(xval-offset,yval,str(int(P_detection[y_id,x_id])),color=col,fontsize=ft)

        ax2 = fig2.add_subplot(3,3,i+1+(3*(id_mel)))
        im = ax2.pcolor(T_x,T_y,I, cmap=cmap, vmin=0, vmax=2)
        cb = fig2.colorbar(im, ax=ax2)
        cb.set_label("$\phi$ (a.u)",fontsize=ft_label)
        ax2.set_xticks(x)  # Set x ticks to vector values
        ax2.set_yticks(y)  # Set y ticks to vector values
        ax2.set_xlabel("Source detector separation (cm)",fontsize=ft)
        ax2.set_ylabel("Melanosome volume fraction (%)",fontsize=ft)




        ax3 = fig3.add_subplot(3,3,i+1+(3*(id_mel)))
        ax3.set_title("Placenta depth "+str(dist_to_plancenta_subjects_mm[i])+" mm",fontsize=ft)
        im = ax3.pcolor(T_x,T_y,P_noise, vmin=0, vmax=100, cmap=cmap)
        cb = fig3.colorbar(im, ax=ax3)
        cb.set_label("Monte Carlo noise probability (%)",fontsize=ft_label)
        ax3.set_xticks(x)  # Set x ticks to vector values
        ax3.set_yticks(y)  # Set y ticks to vector values
        ax3.set_xlabel("Source detector separation (cm)",fontsize=ft)
        ax3.set_ylabel("Placenta blood volume ($\mu$M)",fontsize=ft)
plt.show()


# mup_min=np.array([0.01481771, -0.46458202, -0.18996804, -0.77454701, -0.24943377,-0.86297867])


## Process scanning probability with ultrasound measurements

device_name = "Mini CYRIL"
# device_name = "CYRIL"


# title_distrib = "sensitivity indexes"
title_distrib = "scanning probability"




# Extract thickness from segmentation files
info_thickness = info_subject()
info_thickness = read_thickness_data(info_thickness)
info_thickness = read_thickness_data(info_thickness,col='.1')
info_thickness = read_thickness_data(info_thickness,col='.2')


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



#load tissue sensitivity
data = np.load(main_path+"simulations/tissue_sensitivity.npz")
Tissue_sensitivity_indexes = data['Tissue_sensitivity_indexes']
skin_thickness_array = data['skin_thickness_array']
adipose_thickness_array = data['adipose_thickness_array']
muscle_thickness_array = data['muscle_thickness_array']
SatO2_muscle_array = data['SatO2_muscle_array']
SatO2_placenta_array = data['SatO2_placenta_array']
f_melanosome = data['f_melanosome']
HbT_muscle_array = data['HbT_muscle_array']
HbT_placenta_array = data['HbT_placenta_array']
SD_separation_cm = data['SD_separation_cm']

if device_name == "Mini CYRIL":
    SD_separation_cm = SD_separation_cm[2:5]
    int_time_array = np.array([1,5,10])
else:
    int_time_array = np.array([1])


for int_time_s in int_time_array:
    for binning in np.array([1,10]):
        print("it ",int_time_s,"binning",binning)

        Proba_matrix = Tissue_sensitivity_indexes.copy()

        if device_name == "Mini CYRIL":
            Proba_matrix = Proba_matrix[:,:,:,:,:,:,:,:,2:5,:]


        # Proba_matrix = Scanning_probability.copy()
        if title_distrib == "scanning probability":

            #Load detection proba
            data = np.load(main_path+"simulations/detection_proba/detection_probability_"+device_name+"_TI_"+str(int_time_s)+"s_binning_"+str(binning)+"_camera_model_0.npz")
            Detection_probability = data['Detection_probability'] #convert in %
            detector_signal = data['detector_signal']
            noise_level = data['noise_level']

            for i in range(Proba_matrix.shape[9]):
                Proba_matrix[:,:,:,:,:,:,:,:,:,i]*=Detection_probability


        #Create maps of tissue sensitivity in function of SatO2_muscle and SatO2_placenta
        #knowing skin tones, thicknesses (skin, adipose, muscle)

        #Init maps
        S_skin = np.zeros((SD_separation_cm.shape[0], t_skin.shape[0], SatO2_muscle_array.shape[0], SatO2_placenta_array.shape[0],HbT_muscle_array.shape[0],HbT_placenta_array.shape[0]))
        S_adipose = np.zeros(S_skin.shape)
        S_muscle = np.zeros(S_skin.shape)
        S_placenta = np.zeros(S_skin.shape)
        det_proba = np.zeros(S_skin.shape)

        dr = np.zeros(S_skin.shape)

        #Loop on detectors
        for d in range(SD_separation_cm.shape[0]):
            #Loop on subject ID
            for i in range(t_skin.shape[0]):

                for sat_m in range(SatO2_muscle_array.shape[0]):
                    for sat_p in range(SatO2_placenta_array.shape[0]):
                        for HbT_m in range(HbT_muscle_array.shape[0]):
                            for HbT_p in range(HbT_placenta_array.shape[0]):

                                #Interpolation functions
                                interp_skin = interpolate.RegularGridInterpolator(
                                (skin_thickness_array, adipose_thickness_array,
                                muscle_thickness_array,f_melanosome),
                                Proba_matrix[:,:,:,sat_m,sat_p,:,HbT_m,HbT_p,d,0],
                                bounds_error=False, fill_value=None)

                                interp_adipose = interpolate.RegularGridInterpolator(
                                (skin_thickness_array, adipose_thickness_array,
                                muscle_thickness_array,f_melanosome),
                                Proba_matrix[:,:,:,sat_m,sat_p,:,HbT_m,HbT_p,d,1],
                                bounds_error=False, fill_value=None)

                                interp_muscle = interpolate.RegularGridInterpolator(
                                (skin_thickness_array, adipose_thickness_array,
                                muscle_thickness_array,f_melanosome),
                                Proba_matrix[:,:,:,sat_m,sat_p,:,HbT_m,HbT_p,d,2],
                                bounds_error=False, fill_value=None)

                                interp_placenta = interpolate.RegularGridInterpolator(
                                (skin_thickness_array, adipose_thickness_array,
                                muscle_thickness_array,f_melanosome),
                                Proba_matrix[:,:,:,sat_m,sat_p,:,HbT_m,HbT_p,d,3],
                                bounds_error=False, fill_value=None)


                                S_skin[d,i, sat_m,sat_p,HbT_m,HbT_p] = interp_skin(((t_skin[i]*10),
                                (t_adipose[i]*10), (t_muscle[i]*10), f_mel[i]))

                                S_adipose[d,i, sat_m,sat_p,HbT_m,HbT_p] = interp_adipose(((t_skin[i]*10),
                                (t_adipose[i]*10), (t_muscle[i]*10), f_mel[i]))

                                S_muscle[d,i, sat_m,sat_p,HbT_m,HbT_p] = interp_muscle(((t_skin[i]*10),
                                (t_adipose[i]*10), (t_muscle[i]*10), f_mel[i]))

                                S_placenta[d,i, sat_m,sat_p,HbT_m,HbT_p] = interp_placenta(((t_skin[i]*10),
                                (t_adipose[i]*10), (t_muscle[i]*10), f_mel[i]))


                                if title_distrib == "scanning probability":
                                    interp_dr = interpolate.RegularGridInterpolator(
                                    (skin_thickness_array, adipose_thickness_array,
                                    muscle_thickness_array,f_melanosome),
                                    detector_signal[:,:,:,sat_m,sat_p,:,HbT_m,HbT_p,d],
                                    bounds_error=False, fill_value=None)

                                    dr[d,i, sat_m,sat_p,HbT_m,HbT_p] = interp_dr(((t_skin[i]*10),
                                    (t_adipose[i]*10), (t_muscle[i]*10), f_mel[i]))

                                    interp_det_proba = interpolate.RegularGridInterpolator(
                                    (skin_thickness_array, adipose_thickness_array,
                                    muscle_thickness_array,f_melanosome),
                                    Detection_probability[:,:,:,sat_m,sat_p,:,HbT_m,HbT_p,d],
                                    bounds_error=False, fill_value=None)

                                    det_proba[d,i, sat_m,sat_p,HbT_m,HbT_p] = interp_det_proba(((t_skin[i]*10), (t_adipose[i]*10), (t_muscle[i]*10), f_mel[i]))



        S_skin[S_skin<0] = 0
        S_skin[S_skin>1] = 1
        S_adipose[S_adipose<0] = 0
        S_adipose[S_adipose>1] = 1
        S_muscle[S_muscle<0] = 0
        S_muscle[S_muscle>1] = 1

        S_placenta[S_placenta<0]=0
        S_placenta[S_placenta>1]=1



        if title_distrib == "scanning probability":
            np.savez(main_path+"simulations/scanning_proba/"+device_name+"_ti_"+str(int_time_s)+"_binning_"+str(binning)+".npz",
            S_skin = S_skin,
            S_adipose = S_adipose,
            S_muscle = S_muscle,
            S_placenta = S_placenta,
            t_skin = t_skin,
            t_adipose = t_adipose,
            t_muscle = t_muscle,
            Fitzpatrick = Fitzpatrick,
            ID_subjects = ID_subjects,
            Gestation_subjects = Gestation_subjects,
            SD_separation_cm = SD_separation_cm,
            detector_signal = dr,
            Detection_proba = det_proba,
            noise_level = noise_level)
        else:
            np.savez(main_path+"simulations/scanning_proba/sensitivity_indexes.npz",
            S_skin = S_skin,
            S_adipose = S_adipose,
            S_muscle = S_muscle,
            S_placenta = S_placenta,
            t_skin = t_skin,
            t_adipose = t_adipose,
            t_muscle = t_muscle,
            Fitzpatrick = Fitzpatrick,
            ID_subjects = ID_subjects,
            Gestation_subjects = Gestation_subjects,
            SD_separation_cm = SD_separation_cm)

## Process Mini CYRIL proba with ultrasound measurements v2

device_name = "Mini CYRIL"


# Extract thickness from segmentation files
info_thickness = info_subject()
info_thickness = read_thickness_data(info_thickness)
info_thickness = read_thickness_data(info_thickness,col='.1')
info_thickness = read_thickness_data(info_thickness,col='.2')


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




#remove subject with t_skin>0.5
idx = t_skin<=0.3
t_skin = t_skin[idx]
t_adipose = t_adipose[idx]
t_muscle = t_muscle[idx]
Fitzpatrick = Fitzpatrick[idx]
Gestation_subjects = Gestation_subjects[idx]

#Convert Fitzpatrick scale in f_mel
f_mel = convert_Fitzpatrick_scale_f_mel(Fitzpatrick)

dist_to_placenta = t_skin+t_adipose+t_muscle

#Convert thickness in mm
t_skin *=10
t_adipose *=10
t_muscle *=10
dist_to_placenta *=10

#load tissue sensitivity
data = np.load(main_path+"simulations/tissue_sensitivity.npz")
Tissue_sensitivity_indexes = data['Tissue_sensitivity_indexes']
skin_thickness_array = data['skin_thickness_array']
adipose_thickness_array = data['adipose_thickness_array']
muscle_thickness_array = data['muscle_thickness_array']
SatO2_muscle_array = data['SatO2_muscle_array']
SatO2_placenta_array = data['SatO2_placenta_array']
f_melanosome = data['f_melanosome']
HbT_muscle_array = data['HbT_muscle_array']
HbT_placenta_array = data['HbT_placenta_array']
SD_separation_cm = data['SD_separation_cm']


idx_SD = np.array([2,3,4])
SD_separation_cm = SD_separation_cm[idx_SD]
int_time_array = np.array([1,5,10])
Sensitivity_Placenta = Tissue_sensitivity_indexes[:,:,:,:,:,:,:,:,idx_SD,-1]



for int_time_s in int_time_array:
    for binning in np.array([1,10]):
        print("it ",int_time_s,"binning",binning)

        #Load detection proba
        data = np.load(main_path+"simulations/detection_proba/detection_probability_"+device_name+"_TI_"+str(int_time_s)+"s_binning_"+str(binning)+"_camera_model_0.npz")
        Detection_probability = data['Detection_probability'] #convert in %


        S_placenta = []
        P_detection = []
        P_scanning = []

        #Loop over subjects
        for id_subject in range(t_skin.shape[0]):
            #Get idx subject in LUT
            id_skin = np.argmin(np.abs(skin_thickness_array - t_skin[id_subject]))
            id_adipose = np.argmin(np.abs(adipose_thickness_array - t_adipose[id_subject]))
            id_muscle = np.argmin(np.abs(muscle_thickness_array - t_muscle[id_subject]))
            id_mel = np.argmin(np.abs(f_melanosome-f_mel[id_subject]))

            #get proba
            S_placenta_subject = Sensitivity_Placenta[id_skin,id_adipose,id_muscle, :, :, id_mel, :, :, :]
            P_detection_subject = Detection_probability[id_skin,id_adipose,id_muscle, :, :, id_mel, :, :, :]
            Scanning_proba_subject = np.multiply(S_placenta_subject,P_detection_subject)

            S_placenta.append(S_placenta_subject)
            P_detection.append(P_detection_subject)
            P_scanning.append(Scanning_proba_subject)


        S_placenta = np.asarray(S_placenta)
        P_detection = np.asarray(P_detection)
        P_scanning = np.asarray(P_scanning)


        np.savez(main_path+"simulations/scanning_proba/"+device_name+"_ti_"+str(int_time_s)+"_binning_"+str(binning)+".npz",
        S_placenta = S_placenta,
        P_detection = P_detection,
        P_scanning = P_scanning,
        t_skin = t_skin,
        t_adipose = t_adipose,
        t_muscle = t_muscle,
        f_mel = f_mel,
        Fitzpatrick = Fitzpatrick,
        ID_subjects = ID_subjects,
        Gestation_subjects = Gestation_subjects,
        SD_separation_cm = SD_separation_cm,
        SatO2_muscle_array = SatO2_muscle_array,
        SatO2_placenta_array = SatO2_placenta_array,
        HbT_muscle_array = HbT_muscle_array,
        HbT_placenta_array=HbT_placenta_array)


## Create Sensitivity xlsx file
Proba = np.copy(S_placenta)
Proba = 100*Proba
Proba[Proba<0] = 0
Proba[Proba>100] = 100

HbT_m = 25
HbT_p = 25
id_HbT_m = np.where(HbT_muscle_array == HbT_m)[0].item()
id_HbT_p = np.where(HbT_placenta_array == HbT_p)[0].item()


if device_name == "CYRIL":
    df = pd.DataFrame(
    {
    "ID subjects": ID_subjects,
    "Gestation week": Gestation_subjects,
    "Mean sensitivity placenta (2 cm)" : Proba[1,:,:,:,id_HbT_m,id_HbT_p].mean(axis=(1,2)),
    # "Std sensitivity placenta (2 cm)" : Proba[1,:,:,:,id_HbT_m,id_HbT_p].std(axis=(1,2)),

    "Mean sensitivity placenta (3 cm)" : Proba[2,:,:,:,id_HbT_m,id_HbT_p].mean(axis=(1,2)),
    # "Std sensitivity placenta (3 cm)" : Proba[2,:,:,:,id_HbT_m,id_HbT_p].std(axis=(1,2)),

    "Mean sensitivity placenta (4 cm)" : Proba[3,:,:,:,id_HbT_m,id_HbT_p].mean(axis=(1,2)),
    # "Std sensitivity placenta (4 cm)" : Proba[3,:,:,:,id_HbT_m,id_HbT_p].std(axis=(1,2)),
    }
    )
else:
    df = pd.DataFrame(
    {
    "ID subjects": ID_subjects,
    "Gestation week": Gestation_subjects,

    "Mean sensitivity placenta (3 cm)" : Proba[0,:,:,:,id_HbT_m,id_HbT_p].mean(axis=(1,2)),
    # "Std sensitivity placenta (3 cm)" : Proba[0,:,:,:,id_HbT_m,id_HbT_p].std(axis=(1,2)),

    "Mean sensitivity placenta (4 cm)" : Proba[1,:,:,:,id_HbT_m,id_HbT_p].mean(axis=(1,2)),
    # "Std sensitivity placenta (4 cm)" : Proba[1,:,:,:,id_HbT_m,id_HbT_p].std(axis=(1,2)),

    "Mean sensitivity placenta (5 cm)" : Proba[2,:,:,:,id_HbT_m,id_HbT_p].mean(axis=(1,2)),
    # "Std sensitivity placenta (5 cm)" : Proba[2,:,:,:,id_HbT_m,id_HbT_p].std(axis=(1,2)),
    }
    )

df.to_excel('/home/caredda/temp/Placenta_sensitivity_'+device_name+'_TI'+str(integration_time_s)+'.xlsx')




## Plot proba for all subjects


device_name = "Mini CYRIL"

plt.close('all')
fig = plt.figure(tight_layout=True)
gs = gridspec.GridSpec(3, 1)

ft = 16
c = plt.rcParams['axes.prop_cycle'].by_key()['color']


S_placenta_m = []
P_detection_m = []
P_scanning_m = []
x_ticks = []
x_labels = []
quantiles_vec = []
colors=[]
patches = []

# ref:Non-invasive monitoring of blood oxygenation in human placentas via concurrent diffuse optical spectroscopy and ultrasound imaging
id_HbT_placenta = 2 # 35 µM
id_SatO2_placenta = 2 #80%

id_HbT_muscle = 2 # 35 µM
id_SatO2_muscle = 1 #60% ref: 68% in Changes in Muscle Oxygen Saturation Measured Using Wireless Near-Infrared Spectroscopy in Resistance Training: A Systematic Review





id = 0
for d in np.arange(3):
    for int_time_s in np.array([1,5,10]):
        for binning in np.array([1,10]):



            #Load data
            filename = main_path+"simulations/scanning_proba/"+device_name+"_ti_"+str(int_time_s)+"_binning_"+str(binning)+".npz"
            data = np.load(filename)
            S_placenta = data['S_placenta']
            P_detection = data['P_detection']
            P_scanning = data['P_scanning']
            t_skin = data['t_skin']
            t_adipose = data['t_adipose']
            t_muscle = data['t_muscle']
            f_mel = data['f_mel']
            Fitzpatrick = data['Fitzpatrick']
            ID_subjects = data['ID_subjects']
            Gestation_subjects = data['Gestation_subjects']
            SD_separation_cm = data['SD_separation_cm'].astype(np.uint8)
            SatO2_muscle_array = data['SatO2_muscle_array']
            SatO2_placenta_array = data['SatO2_placenta_array']
            HbT_muscle_array = data['HbT_muscle_array']
            HbT_placenta_array= data['HbT_placenta_array']

            S_placenta_m.append(S_placenta[:,id_SatO2_muscle,id_SatO2_placenta,id_HbT_muscle,id_HbT_placenta, d])
            P_detection_m.append(P_detection[:,id_SatO2_muscle,id_SatO2_placenta,id_HbT_muscle,id_HbT_placenta, d])
            P_scanning_m.append(P_scanning[:,id_SatO2_muscle,id_SatO2_placenta,id_HbT_muscle,id_HbT_placenta, d])
            # S_placenta_m.append(S_placenta.mean(axis=(1,2,3,4))[:,d])
            # P_detection_m.append(P_detection.mean(axis=(1,2,3,4))[:,d])
            # P_scanning_m.append(P_scanning.mean(axis=(1,2,3,4))[:,d])


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

#Placenta sensitivity
S_placenta_m = np.asarray(S_placenta_m).T
S_placenta_m = S_placenta_m[:,np.array([0,6,12])]

ax = fig.add_subplot(gs[0,0])
ax.set_title("Placenta sensitivity",fontsize = ft)
plots = ax.violinplot(S_placenta_m*100, showmedians=True, showmeans=False,quantiles=[[0.25,0.75],[0.25,0.75],[0.25,0.75]])
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
P_detection_m = np.asarray(P_detection_m).T
ax = fig.add_subplot(gs[1,0])
ax.set_title("Detection probability",fontsize = ft)
plots = ax.violinplot(P_detection_m*100, showmedians=True, showmeans=False,quantiles=quantiles_vec)
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)
plots['cmedians'].set_colors(c[id])
plots['cquantiles'].set_color('r')
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
P_scanning_m = np.asarray(P_scanning_m).T
ax = fig.add_subplot(gs[2,0])
ax.set_title("Placenta scanning probability",fontsize = ft)
plots = ax.violinplot(P_scanning_m*100, showmedians=True, showmeans=False,quantiles=quantiles_vec)
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)
plots['cmedians'].set_colors(c[id])
plots['cquantiles'].set_color('r')
ax.set_xticks(np.array([3.5,9.5,15.5]))
ax.set_xticklabels(np.array(["3","4","5"]),fontsize=ft)
ax.set_ylabel("Probability (%)",fontsize = ft)
ax.set_xlabel("Source-detector separation (cm)",fontsize = ft)
ax.axvline(6.5,0,100,color='k',linewidth=3)
ax.axvline(12.5,0,100,color='k',linewidth=3)
ax.grid()
ax.set_ylim(0,100)
# put those patched as legend-handles into the legend
ax.legend(handles=patches,loc="best", borderaxespad=0. ,fontsize = ft)


        # id += 1
plt.show()



## Plot proba for all subjects - ECBO Plot


device_name = "Mini CYRIL"

plt.close('all')
fig = plt.figure(tight_layout=True)
gs = gridspec.GridSpec(2, 1)

ft = 16
c = plt.rcParams['axes.prop_cycle'].by_key()['color']
plt.rcParams['font.size'] = ft
S_placenta_m = []
P_detection_m = []
P_scanning_m = []
x_ticks = []
x_labels = []
quantiles_vec = []
colors=[]
patches = []

# ref:Non-invasive monitoring of blood oxygenation in human placentas via concurrent diffuse optical spectroscopy and ultrasound imaging
id_HbT_placenta = 2 # 35 µM
id_SatO2_placenta = 2 #80%

id_HbT_muscle = 2 # 35 µM
id_SatO2_muscle = 1 #60% ref: 68% in Changes in Muscle Oxygen Saturation Measured Using Wireless Near-Infrared Spectroscopy in Resistance Training: A Systematic Review


int_time_array = np.array([5])
binning_array = np.array([1,10])


id = 0
for d in np.arange(3):
    for int_time_s in int_time_array:
        for binning in binning_array:



            #Load data
            filename = main_path+"simulations/scanning_proba/"+device_name+"_ti_"+str(int_time_s)+"_binning_"+str(binning)+".npz"
            data = np.load(filename)
            S_placenta = data['S_placenta']
            P_detection = data['P_detection']
            P_scanning = data['P_scanning']
            t_skin = data['t_skin']
            t_adipose = data['t_adipose']
            t_muscle = data['t_muscle']
            f_mel = data['f_mel']
            Fitzpatrick = data['Fitzpatrick']
            ID_subjects = data['ID_subjects']
            Gestation_subjects = data['Gestation_subjects']
            SD_separation_cm = data['SD_separation_cm'].astype(np.uint8)
            SatO2_muscle_array = data['SatO2_muscle_array']
            SatO2_placenta_array = data['SatO2_placenta_array']
            HbT_muscle_array = data['HbT_muscle_array']
            HbT_placenta_array= data['HbT_placenta_array']

            S_placenta_m.append(S_placenta[:,id_SatO2_muscle,id_SatO2_placenta,id_HbT_muscle,id_HbT_placenta, d])
            P_detection_m.append(P_detection[:,id_SatO2_muscle,id_SatO2_placenta,id_HbT_muscle,id_HbT_placenta, d])
            P_scanning_m.append(P_scanning[:,id_SatO2_muscle,id_SatO2_placenta,id_HbT_muscle,id_HbT_placenta, d])
            # S_placenta_m.append(S_placenta.mean(axis=(1,2,3,4))[:,d])
            # P_detection_m.append(P_detection.mean(axis=(1,2,3,4))[:,d])
            # P_scanning_m.append(P_scanning.mean(axis=(1,2,3,4))[:,d])


for i in range(len(P_scanning_m)):
    quantiles_vec.append([0.25,0.75])

for i in range(3):

    colors.append(c[0])
    colors.append(c[1])


for int_time_s in int_time_array:
    id = 0
    for binning in binning_array:
        patches.append(mpatches.Patch(color=c[id], label="Integration time "+str(int_time_s)+" s, binning "+str(binning)))
        id += 1

#Placenta sensitivity
S_placenta_m = np.asarray(S_placenta_m).T
S_placenta_m = S_placenta_m[:,np.array([0,2,4])]

ax = fig.add_subplot(gs[0,0])
ax.set_title("Placenta sensitivity",fontsize = ft)
plots = ax.violinplot(S_placenta_m*100, showmedians=True, showmeans=False,quantiles=[[0.25,0.75],[0.25,0.75],[0.25,0.75]])
plots['cmedians'].set_colors('k')
plots['cquantiles'].set_color('r')
ax.set_xticks(np.array([1,2,3]))
ax.set_xticklabels(np.array(["3","4","5"]),fontsize=ft)
ax.set_ylabel("Probability (%)",fontsize = ft)
ax.set_xlabel("Source-detector separation (cm)",fontsize = ft)

# for i in range(len(labels)-1):
#     ax.axvline(SD_cm.shape[0]+0.5+SD_cm.shape[0]*i)
ax.grid()
ax.set_ylim(0,100)


# #Detection probability
# P_detection_m = np.asarray(P_detection_m).T
# ax = fig.add_subplot(gs[1,0])
# ax.set_title("Detection probability",fontsize = ft)
# plots = ax.violinplot(P_detection_m*100, showmedians=True, showmeans=False,quantiles=quantiles_vec)
# for pc, color in zip(plots['bodies'], colors):
#     pc.set_facecolor(color)
# plots['cmedians'].set_colors(c[id])
# plots['cquantiles'].set_color('r')
# ax.set_xticks(np.array([3.5,9.5,15.5]))
# ax.set_xticklabels(np.array(["3","4","5"]),fontsize=ft)
# ax.set_ylabel("Probability (%)",fontsize = ft)
# ax.set_xlabel("Source-detector separation (cm)",fontsize = ft)
# ax.axvline(6.5,0,100,color='k',linewidth=3)
# ax.axvline(12.5,0,100,color='k',linewidth=3)
# ax.grid()
# ax.set_ylim(0,100)
# ax.legend(handles=patches,loc="best", borderaxespad=0. ,fontsize = ft)


#Scanning probability
P_scanning_m = np.asarray(P_scanning_m).T
ax = fig.add_subplot(gs[1,0])
ax.set_title("Placenta scanning probability",fontsize = ft)
plots = ax.violinplot(P_scanning_m*100, showmedians=True, showmeans=False,quantiles=quantiles_vec)
for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)
plots['cmedians'].set_colors('k')
plots['cquantiles'].set_color('r')
ax.set_xticks(np.array([1.5,3.5,5.5]))
ax.set_xticklabels(np.array(["3","4","5"]),fontsize=ft)
ax.set_ylabel("Probability (%)",fontsize = ft)
ax.set_xlabel("Source-detector separation (cm)",fontsize = ft)
ax.axvline(2.5,0,100,color='k',linewidth=3)
ax.axvline(4.5,0,100,color='k',linewidth=3)
ax.grid()
ax.set_ylim(0,100)
# put those patched as legend-handles into the legend
ax.legend(handles=patches,loc="best", borderaxespad=0. ,fontsize = ft)


        # id += 1
plt.show()

## Study placenta scanning probability as a function of integration time and temporal binning [v2]


# Extract thickness from segmentation files
info_thickness = info_subject()
info_thickness = read_thickness_data(info_thickness)
info_thickness = read_thickness_data(info_thickness,col='.1')
info_thickness = read_thickness_data(info_thickness,col='.2')


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

id_nan = ~ np.isnan(t_muscle)
t_skin = t_skin[id_nan]
t_adipose = t_adipose[id_nan]
t_muscle = t_muscle[id_nan]




ft_legend = 14
ft = 16
c = plt.rcParams['axes.prop_cycle'].by_key()['color']

plt.close('all')

##

percentile = 25
val_skin_mm = np.percentile(t_skin,percentile)*10
val_adipose_mm = np.percentile(t_adipose,percentile)*10
val_muscle_mm = np.percentile(t_muscle,percentile)*10

print("skin",val_skin_mm/10)
print("Adipose",val_adipose_mm/10)
print("Muscle",val_muscle_mm/10)


dist_to_muscle = val_skin_mm + val_adipose_mm
dist_to_placenta = val_skin_mm + val_adipose_mm + val_muscle_mm



Scanning_CYRIL = []
Scanning_Mini_CYRIL = []

Scanning_CYRIL_muscle = []
Scanning_Mini_CYRIL_muscle  = []
labels = []


for device_name in np.array(["CYRIL", "Mini CYRIL"]):

    if device_name == "CYRIL":
        int_time_array = np.array([1])

    if device_name == "Mini CYRIL":
        int_time_array = np.array([1,5,10])

    for int_time_s in int_time_array:
        for binning in np.array([1,10]):
            _S_muscle, _S_placenta = calculate_scanning_proba(val_skin_mm, val_adipose_mm, val_muscle_mm, device_name,binning, int_time_s)


            if device_name == "CYRIL":
                Scanning_CYRIL.append(100*_S_placenta)
                Scanning_CYRIL_muscle.append(100*_S_muscle)

            if device_name == "Mini CYRIL":
                Scanning_Mini_CYRIL.append(100*_S_placenta)
                Scanning_Mini_CYRIL_muscle.append(100*_S_muscle)

            labels.append(device_name+" (IT = "+str(int_time_s)+" s, temporal binning "+str(binning)+")")














#Rearange data
Scanning_proba = []
Scanning_proba_muscle = []
colors = []
quantiles = []

for i in range(2):
    for j in range(len(Scanning_CYRIL)):
        Scanning_proba.append(Scanning_CYRIL[j][i,:,:,:,:].ravel())
        Scanning_proba_muscle.append(Scanning_CYRIL_muscle[j][i,:,:,:,:].ravel())
        quantiles.append(np.array([0.25,0.75]))
        colors.append(c[j])
    for j in range(len(Scanning_Mini_CYRIL)):
        Scanning_proba.append(np.zeros((2,)))
        Scanning_proba_muscle.append(np.zeros((2,)))
        quantiles.append(np.array([0.25,0.75]))
        colors.append(c[j+2])

for i in range(3):
    for j in range(len(Scanning_CYRIL)):
        Scanning_proba.append(Scanning_CYRIL[j][2+i,:,:,:,:].ravel())
        Scanning_proba_muscle.append(Scanning_CYRIL_muscle[j][2+i,:,:,:,:].ravel())
        quantiles.append(np.array([0.25,0.75]))
        colors.append(c[j])
    for j in range(len(Scanning_Mini_CYRIL)):
        Scanning_proba.append(Scanning_Mini_CYRIL[j][i,:,:,:,:].ravel())
        Scanning_proba_muscle.append(Scanning_Mini_CYRIL_muscle[j][i,:,:,:,:].ravel())

        quantiles.append(np.array([0.25,0.75]))
        colors.append(c[j+2])

# Plots
plt.close('all')





ft_legend2 = 12



plt.figure()

#Muscle scanning proba
ax2 = plt.subplot(212)
ax2.set_title("Muscle scanning probability (Distance to muscle: "+str(int(dist_to_muscle*100)/100)+" mm)",fontsize=ft)

plots = ax2.violinplot(Scanning_proba_muscle,showmedians=True, showmeans=False,quantiles=quantiles)


for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)
plots['cmedians'].set_colors(colors)
plots['cquantiles'].set_color('y')


for i in range(4):
    plt.axvline(len(labels)+0.5+len(labels)*i)

# create a patch (proxy artist) for every color
patches = []
for i in range(len(labels)):
    patches.append(mpatches.Patch(color=colors[i], label=labels[i]))
plt.legend(handles=patches,loc="upper left", borderaxespad=0.,fontsize = ft_legend2 )

ax2.grid()
ax2.set_xlabel("Source detector separation (mm)",fontsize=ft)
ax2.set_ylabel("Probability (%)",fontsize=ft)
ax2.set_xticks(np.arange(len(labels)/2+0.5,5*(len(labels)+0.5),len(labels)),np.array(["10","20","30","40","50"]))
# ax2.set_ylim(ax2.get_ylim())







ax1 = plt.subplot(211)
ax1.set_title("Placenta scanning probability (Distance to placenta: "+str(int(dist_to_placenta*100)/100)+" mm)",fontsize=ft)

plots = ax1.violinplot(Scanning_proba,showmedians=True, showmeans=False,quantiles=quantiles)


for pc, color in zip(plots['bodies'], colors):
    pc.set_facecolor(color)
plots['cmedians'].set_colors(colors)
plots['cquantiles'].set_color('y')


for i in range(4):
    plt.axvline(len(labels)+0.5+len(labels)*i)

# create a patch (proxy artist) for every color
patches = []
for i in range(len(labels)):
    patches.append(mpatches.Patch(color=colors[i], label=labels[i]))
plt.legend(handles=patches,loc="upper left", borderaxespad=0.,fontsize = ft_legend2)

ax1.grid()
ax1.set_xlabel("Source detector separation (mm)",fontsize=ft)
ax1.set_ylabel("Probability (%)",fontsize=ft)
ax1.set_xticks(np.arange(len(labels)/2+0.5,5*(len(labels)+0.5),len(labels)),np.array(["10","20","30","40","50"]))
# ax1.set_ylim(ax2.get_ylim())






plt.show()



## Plot sensitivity profiles


path = main_path + "simulations/output_lookup_table_skin/"

#pos layers
Thickness_skin = 1 #1,2,3,5
Thickness_adipose = 4
Thickness_muscle = 10




ft_txt = 15
lw = 3
plt.close('all')
plt.figure()

plt.rcParams.update({'font.size': ft_txt})

data = scipy.io.loadmat(path+"out_St_muscle_0.6_St_placenta_0.8_Thick_skin_"+str(Thickness_skin)+"_Thick_adipose_4_Thick_muscle_10f_mel0.305_HbT_muscle_umol_25_HbT_placenta_umol_25.mat")

sensitivity_map = data['Sensitivity_maps']
#Normalize by te sum for each detector
for d in range(sensitivity_map.shape[0]):
    sensitivity_map[d,:,:,:] /= sensitivity_map[d,:,:,:].sum()

#Convert in percentage
sensitivity_map *= 100



# Source pos
volume_square_size = 200
reso = 1
src_pos = np.array([int(volume_square_size/2 -1),int((volume_square_size - 80/reso)/2 -1),0])


#Define multiple detectors  (values in accordance with VOlume size, se below)
det_pos = np.zeros((SD_separation_cm.shape[0],2),dtype=int)
for i in range(SD_separation_cm.shape[0]):
    det_pos[i,0] = src_pos[0]
    det_pos[i,1] = src_pos[1]+(SD_separation_cm[i])*int(10/reso)


# only use the detectors det_pos
sensitivity_map = sensitivity_map[(SD_separation_cm-1).astype(int),:,:,:]




for id_det,det in enumerate(SD_separation_cm):


    map = sensitivity_map[id_det,det_pos[0,0],src_pos[1]-10:det_pos[id_det,1]+10,1:22].copy()
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


    p = plt.subplot(1,SD_separation_cm.shape[0],id_det+1)
    plt.title("Source-detector separation "+str(int(SD_separation_cm[id_det]))+" cm")

    im = plt.contourf(x_interp, y_interp,map_interp,levels = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100], locator=ticker.LogLocator(),cmap='plasma')

    plt.plot(x_interp,Thickness_skin*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)
    plt.plot(x_interp,Thickness_adipose*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)
    plt.plot(x_interp,Thickness_muscle*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)

    plt.plot(10,0,'ko',markersize=12,label="Source")
    plt.plot(10+det_pos[id_det,1]-src_pos[1],0,'go',markersize=12,label="Detector")

    plt.text(1,Thickness_skin-0.2,"Skin",fontsize=ft_txt)
    plt.text(1,Thickness_adipose-0.2,"Adipose tissue",fontsize=ft_txt)
    plt.text(1,Thickness_muscle-0.2,"Muscle",fontsize=ft_txt)
    plt.text(1,Thickness_muscle+0.5,"Placenta",fontsize=ft_txt)





    c = plt.colorbar(im)
    c.set_label("Sensitivity probability (%)")
    p.invert_yaxis()
    # plt.xticks(x)
    # plt.yticks(y)
    plt.xlabel("Tissue width (mm)")
    plt.ylabel("Tissue depth (mm)")
    plt.legend(loc="best")
plt.show()


plt.figure()

plt.title("Source-detector separation "+str(int(SD_separation_cm[id_det]))+" cm")

im = plt.contourf(x_interp, y_interp,map_interp,levels = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100], locator=ticker.LogLocator(),cmap='plasma')

plt.plot(x_interp,Thickness_skin*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)
plt.plot(x_interp,Thickness_adipose*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)
plt.plot(x_interp,Thickness_muscle*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)

plt.plot(10,0,'ko',markersize=12,label="Source")
plt.plot(10+det_pos[id_det,1]-src_pos[1],0,'go',markersize=12,label="Detector")

plt.text(1,Thickness_skin-0.2,"Skin",fontsize=ft_txt)
plt.text(1,Thickness_adipose-0.2,"Adipose tissue",fontsize=ft_txt)
plt.text(1,Thickness_muscle-0.2,"Muscle",fontsize=ft_txt)
plt.text(1,Thickness_muscle+0.5,"Placenta",fontsize=ft_txt)





c = plt.colorbar(im)
c.set_label("Sensitivity probability (%)")
p.invert_yaxis()
# plt.xticks(x)
# plt.yticks(y)
plt.xlabel("Tissue width (mm)")
plt.ylabel("Tissue depth (mm)")
plt.legend(loc="best")
plt.show()


## Study skin effect

path = main_path + "simulations/output_lookup_table_skin/"

#pos layers
Thickness_skin = 1 #1,2,3,5
Thickness_adipose = 4
Thickness_muscle = 10

f_mel_array = np.array([0.0255,0.305])
SD_separation_cm = np.arange(1,9)


ft = 20
ft_txt = 15
lw = 3
plt.close('all')
plt.figure()

vec_map_interp = []

for f_mel in range(f_mel_array.shape[0]):

    data = scipy.io.loadmat(path+"out_St_muscle_0.6_St_placenta_0.8_Thick_skin_"+str(Thickness_skin)+"_Thick_adipose_4_Thick_muscle_10f_mel"+str(f_mel_array[f_mel])+"_HbT_muscle_umol_25_HbT_placenta_umol_25.mat")

    sensitivity_map = data['Sensitivity_maps']
    #Normalize by te sum for each detector
    for d in range(sensitivity_map.shape[0]):
        sensitivity_map[d,:,:,:] /= sensitivity_map[d,:,:,:].sum()

    #Convert in percentage
    sensitivity_map *= 100

    # Source pos
    volume_square_size = 200
    reso = 1
    src_pos = np.array([int(volume_square_size/2 -1),int((volume_square_size - 80/reso)/2 -1),0])


    #Define multiple detectors  (values in accordance with VOlume size, se below)
    det_pos = np.zeros((SD_separation_cm.shape[0],2),dtype=int)
    for i in range(SD_separation_cm.shape[0]):
        det_pos[i,0] = src_pos[0]
        det_pos[i,1] = src_pos[1]+(i+1)*int(10/reso)






    id_det = 2 # 30 mm
    # map = sensitivity_map[id_det,det_pos[i,0],det_pos[i,1]-50:det_pos[i,1]+50,0:50]
    map = sensitivity_map[id_det,det_pos[0,0],src_pos[1]-10:det_pos[id_det,1]+10,1:22].copy()
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
    p = plt.subplot(1,3,f_mel+1)
    plt.title("Melanosome volume fraction "+str(f_mel_array[f_mel]*100)+"%",fontsize=ft)
    im = plt.contourf(x_interp, y_interp,map_interp,levels = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100], locator=ticker.LogLocator(),cmap='plasma')

    plt.plot(x_interp,Thickness_skin*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)
    plt.plot(x_interp,Thickness_adipose*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)
    plt.plot(x_interp,Thickness_muscle*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)

    plt.plot(10,0,'ko',markersize=12,label="Source")
    plt.plot(40,0,'go',markersize=12,label="Detector")

    plt.text(1,Thickness_skin-0.2,"Skin",fontsize=ft_txt)
    plt.text(1,Thickness_adipose-0.2,"Adipose tissue",fontsize=ft_txt)
    plt.text(1,Thickness_muscle-0.2,"Muscle",fontsize=ft_txt)
    plt.text(1,Thickness_muscle+0.5,"Placenta",fontsize=ft_txt)





    c = plt.colorbar(im)
    c.set_label("Sensitivity probability (%)",fontsize=ft)
    p.invert_yaxis()
    # plt.xticks(x)
    # plt.yticks(y)
    plt.xlabel("Tissue width (mm)",fontsize=ft)
    plt.ylabel("Tissue depth (mm)",fontsize=ft)
    plt.legend(loc="best",fontsize=ft)




#Plot difference map

p = plt.subplot(133)
plt.title("Difference map",fontsize=ft)
im = plt.contourf(x_interp, y_interp,vec_map_interp[1] - vec_map_interp[0])#,levels = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100], locator=ticker.LogLocator(),cmap='plasma')

plt.plot(x_interp,Thickness_skin*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)
plt.plot(x_interp,Thickness_adipose*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)
plt.plot(x_interp,Thickness_muscle*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)

plt.plot(10,0,'ko',markersize=12,label="Source")
plt.plot(40,0,'go',markersize=12,label="Detector")

plt.text(1,Thickness_skin-0.2,"Skin",fontsize=ft_txt)
plt.text(1,Thickness_adipose-0.2,"Adipose tissue",fontsize=ft_txt)
plt.text(1,Thickness_muscle-0.2,"Muscle",fontsize=ft_txt)
plt.text(1,Thickness_muscle+0.5,"Placenta",fontsize=ft_txt)

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
plt.figure()
plt.title("Sensitivity matrix, SD seperation 3 cm",fontsize=ft)
p = plt.subplot(111)
im = plt.contourf(x_interp, y_interp,map_interp,levels = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100], locator=ticker.LogLocator(),cmap='plasma')

plt.plot(x_interp,Thickness_skin*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)
plt.plot(x_interp,Thickness_adipose*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)
plt.plot(x_interp,Thickness_muscle*np.ones(x_interp.shape),'k',linestyle=':',linewidth = lw)

plt.plot(10,0,'ko',markersize=12,label="Source")
plt.plot(40,0,'go',markersize=12,label="Detector")

plt.text(1,Thickness_skin-0.2,"Skin",fontsize=ft_txt)
plt.text(1,Thickness_adipose-0.2,"Adipose tissue",fontsize=ft_txt)
plt.text(1,Thickness_muscle-0.2,"Muscle",fontsize=ft_txt)
plt.text(1,Thickness_muscle+0.5,"Placenta",fontsize=ft_txt)





c = plt.colorbar(im)
c.set_label("Sensitivity probability (%)",fontsize=ft)
p.invert_yaxis()
# plt.xticks(x)
# plt.yticks(y)
plt.xlabel("Tissue width (mm)",fontsize=ft)
plt.ylabel("Tissue depth (mm)",fontsize=ft)
plt.legend(loc="best",fontsize=ft)
plt.show()



## Plot mua placenta muscle

w = np.arange(780,901)
mua_placenta_80 = calculate_mua(w,0.85,0,50e-6,80)
mua_placenta_40 = calculate_mua(w,0.85,0,50e-6,40)

mua_muscle_80 = calculate_mua(w,0.76,0,50e-6,80)
mua_muscle_40 = calculate_mua(w,0.76,0,50e-6,40)

lw = 3

plt.close('all')
plt.figure()
plt.plot(w, mua_placenta_80,linewidth=lw,label = "$\mu_a$ placenta ($SatO_2=80 \%$)")
plt.plot(w, mua_placenta_40,linewidth=lw,label = "$\mu_a$ placenta ($SatO_2=40 \%$)")

# plt.plot(w, mua_muscle_80, 'r:')
# plt.plot(w, mua_muscle_40, 'b:')
# plt.yscale('log')
plt.grid()
plt.xlabel("Wavelength (nm)",fontsize=ft)
plt.ylabel("Absorption coefficient (cm$^{-1}$)",fontsize=ft)
plt.legend(loc="best",fontsize=ft)
plt.title("Absorption coefficient of the placenta for different $SatO_2$ values",fontsize=ft)
plt.xlim(780,900)
plt.show()

## Plot MC noise

files = glob.glob(main_path+"simulations/MC_noise_1GPU/*.mat")
tissue = np.array(["Skin", "Adipose tissue", "Muscle", "Placenta"])
dr_array = []
S_array = []
for f in files:
    data = scipy.io.loadmat(f)
    dr_array.append(data['Diffuse_reflectance'])
    S_array.append(data['Sensitivity_index'])

dr_array = np.squeeze(np.asarray(dr_array))
S_array = np.squeeze(np.asarray(S_array))

ft = 18
width = 2
plt.rcParams.update({'font.size': ft})

#define SD separations in mm
x = np.arange(10,90,10)



#Only select SD separatio between 3 and 5 cm
id_SD = np.array([2,3,4])
x = x[id_SD]
dr_array = dr_array[:,id_SD]
S_array = S_array[:,id_SD,:]





plt.close('all')
plt.subplot(121)
plt.title("Monte Carlo noise on tissue sensitivity",fontsize=ft)
offset=-2
for i in range(4):
    offset+=2
    plt.bar(x+offset,100*S_array.std(axis=0)[:,i],width=width,label=tissue[i])

plt.xlabel("Source detector separation (mm)",fontsize=ft)
plt.ylabel("Monte Carlo noise (%)",fontsize=ft)
plt.xticks(x+2.75,x)
plt.grid()
plt.legend(loc="upper center",fontsize=ft)

plt.subplot(122)
width = 0.2
plt.title("Monte Carlo noise on diffuse reflectance",fontsize=ft)
plt.bar(np.arange(dr_array.shape[1]),dr_array.std(axis=0),width=width)
plt.xlabel("Source detector separation (mm)",fontsize=ft)
plt.ylabel("Monte Carlo noise (a. u.)",fontsize=ft)
plt.xticks(np.arange(dr_array.shape[1]),np.array(["30","40","50"]))
plt.grid()

# plt.yscale("symlog")

plt.show()
