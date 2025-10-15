# process_optical_properties_skin_Fat_muscle_placenta.m
Process optical properties of maternal abdomen (mua, mus, n, g)

# get_epsilon.m
Read absorption and molar extinction coefficient from txt files in the folder "spectra"

# convert_Fitzpatrick_scale_f_mel.m
Convert Fitzpatrick scale (skin tone indicators) into melanosome volume fraction (in %)

# create_meshed_volume_4layers.m
Create a mesh volume having 4 layers

# get_sensitivity_index.m
Calculate tissue senstivity indexes for the 4 tissue layers created with create_meshed_volume_4layers.m

# get_sensitivity_profiles.m
Calculate sensitivity profiles (or Jacobians) in the 4 layers mesh volume
