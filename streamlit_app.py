import streamlit as st
import pandas as pd
import math
import numpy as np
from statistics import median
import plotly.express as px



# Set the title and favicon that appear in the Browser's tab bar.
st.set_page_config(
    page_title='Chlormate: Chlorate Residual Modeller',
    page_icon=':potable_water:', # This is an emoji shortcode. Could be a URL too.
)

# -----------------------------------------------------------------------------
# Draw the actual page

# Set the title that appears at the top of the page.
'''
# :potable_water: Chlormate: Chlorate Modeller

Interactive modeller to inspect variables affecting chlorate formation including temperature, storage time, and hypochlorite concentration.
'''
# provide sliders for selection
min_value = 0
max_value = 151
from_days, to_days = st.slider(
    'Storage Days',
    min_value=min_value,
    max_value=max_value,
    value=[min_value, max_value])
selected_days = np.arange(from_days,to_days,5)

temperatures = np.arange(0,41,5)

parameters = {"Concentration of bulk sodium hypochlorite (wt% Cl2)" : np.arange(3,21.1,1),
"Applied dose of sodium hypochlorite (mg/L)" : np.arange(0.5,7.51,0.5),
"Chlorate content at delivery (wt% Cl2)" : np.arange(2,5.61,0.2),
"pH of bulk chemical" : np.arange(10,14.1,0.5),
"Measured chlorate content (g/L)" : np.arange(0,6.1,0.2)
}

constants = {"Enthalpy": 102.2,
             "Entropy": -55.2,
             "Molar mass of Cl2":70.906,
             "Molar mass of ClO3":83.45,
             "Density of water":1.00,
             "Universal gas constant":8.3145,
             "Volume of an ideal gas":22.41,
             "Intercept" : 0.99881321,
             "Slope" : 0.01652278,
             "Chlorate PCV" : np.array(250),
             "Chlorate TLV" :np.array(150)
}

conversions = {"L to mL": 1000,
               "C to K" : 273.15,
               "J to kJ" :1000,
               "s to days" : 86400
}

if not len(temperatures):
    st.warning("Select at least one tempearture")

selected_temperatures = st.pills(
    f'Average Temperature',
    temperatures,
    selection_mode="multi",
    default = [5,15,25])

temp_arr = np.array(selected_temperatures)

selected_values = {}
for i in parameters:
    if i == 'Measured chlorate content (g/L)':
        selected_values[i] = st.sidebar.select_slider(
            f"{i}" ,
            options = np.round(parameters[i],1),
            value = 0
        )
    else:
        selected_values[i] = st.sidebar.select_slider(
            f"{i}" ,
            options = np.round(parameters[i],1),
            value = median(parameters[i])
        )
    
    
''

# Process the data
specific_gravity = constants['Slope'] * selected_values["Concentration of bulk sodium hypochlorite (wt% Cl2)"] + constants["Intercept"]

free_available_chlorine = (selected_values['Concentration of bulk sodium hypochlorite (wt% Cl2)']/100)*specific_gravity*constants["Density of water"]*conversions["L to mL"]

chlorine_content = free_available_chlorine / constants['Molar mass of Cl2']

hypo_content = chlorine_content

if selected_values["Measured chlorate content (g/L)"] > 0 :
    chlorate_content_g = selected_values["Measured chlorate content (g/L)"]
else:
    chlorate_content_g = selected_values["Chlorate content at delivery (wt% Cl2)"]/100*free_available_chlorine

chlorate_content_mol = chlorate_content_g/constants['Molar mass of ClO3']


# """ make dictionaries for each variable that varies by temperature """
NaOCl_decomposition_rate_constant = {}
NaOCl_ionic_strength = {}
NaOCl_loss_per_second = {}
NaOCl_loss_per_day = {}
Chlorate_formation_per_day = {}
Oxygen_formation_per_day = {}

for i in selected_temperatures:
    j = str(i)
    i = i + conversions['C to K']
    NaOCl_decomposition_rate_constant[j] = (2.083e10 * i * 
                                     math.exp(-((constants['Enthalpy']*conversions['J to kJ'])/(constants['Universal gas constant']*i)))*
                                     math.exp((constants['Entropy']/(constants['Universal gas constant']))))
    NaOCl_ionic_strength[j] = chlorine_content + hypo_content + 10**(-(14-selected_values['pH of bulk chemical']))+chlorate_content_mol
    NaOCl_loss_per_second[j] = 10**((0.143*NaOCl_ionic_strength[j])+np.log10(NaOCl_decomposition_rate_constant[j]))
    NaOCl_loss_per_day[j] = NaOCl_loss_per_second[j] * conversions['s to days']
    Chlorate_formation_per_day[j] = NaOCl_loss_per_day[j]/3.24
    Oxygen_formation_per_day[j] = Chlorate_formation_per_day[j]/8.3

# """ now need to calculate the residual chlorates for each day"""
big_df = pd.DataFrame(columns = ['Time','Temperature','NaOCl','Chlorate','Oxygen','Residual Chlorate (mg/L)','Residual Chlorate (ug/L)'])
for i in selected_temperatures:
    j = str(i)
    df = pd.DataFrame(selected_days, columns = ['Time'])
    df['Temperature'] = i
    df['NaOCl'] = ((hypo_content/(NaOCl_loss_per_day[j]*hypo_content*df['Time']+1))*
                   (constants['Molar mass of Cl2']/(specific_gravity*10)))
    df['Chlorate'] = ((chlorate_content_mol+Chlorate_formation_per_day[j]*
                      ((hypo_content/(NaOCl_loss_per_day[j]*hypo_content*df['Time']+1)))**2*df['Time'])*
                      constants['Molar mass of ClO3'])
    df['Oxygen'] = ((Oxygen_formation_per_day[j]*(hypo_content/(NaOCl_loss_per_day[j]*hypo_content*df['Time']+1))**2*df['Time'])
                     *constants['Volume of an ideal gas'])
    df['Residual Chlorate (mg/L)'] = (df['Chlorate']/(df['NaOCl']*10*specific_gravity)*selected_values['Applied dose of sodium hypochlorite (mg/L)']
    )
    df['Residual Chlorate (ug/L)'] = df['Residual Chlorate (mg/L)'] * 1000    
    big_df = pd.concat([big_df,df])

st.header(('Residual Chlorate Levels'), divider='gray')

# st.line_chart(
#     big_df,
#     x='Time',
#     y='Residual Chlorate (ug/L)',
#     color='Temperature',

# )

# fig,ax = plt.subplots()
# for i in selected_temperatures:
#     selected_df = df[df['Temperature']==int(i)]
#     ax.plot(selected_df['Time'],selected_df['Residual Chlorate (ug/L)'], c = 'red')
# st.pyplot(fig)


fig = px.line(
    big_df,
    x="Time",
    y="Residual Chlorate (ug/L)",
    color="Temperature",
    markers = True,
    color_discrete_sequence=px.colors.qualitative.Set2
)

fig.add_hline(y=250, line_width=3, line_dash="dash", line_color="red",layer = 'below')
st.plotly_chart(fig, use_container_width=True)



st.dataframe(big_df.sort_values(by = ['Time','Temperature'],ascending=True))