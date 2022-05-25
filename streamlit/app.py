from io import StringIO
import streamlit as st
from PIL import Image
import pickle
import base64
import pandas as pd
import numpy as np


# File download
def filedownload(df):
    csv = df.to_csv()
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

def taxa_extractor(input_seq,level=3):
    input_df=pd.read_csv("Profile_sequences_5.csv")
    input_taxa=pd.read_csv("Taxa_info.csv")
    from scipy.spatial import distance

    distance_holder=[]
    for i in range(input_df.shape[0]):
        distance_holder.append(distance.euclidean(input_seq, input_df.iloc[i,:].values[:-1]))
    return(input_taxa.iloc[np.argmin(distance_holder),level])


def significance(input_value):
    input_df=pd.read_csv("baseline_PKJK.csv")
    return(int(len(np.where(input_value > input_df.Value)[0])/len(input_df.Value)*100))

def processor(file_input, label="labels"):
    labels=[]
    sequences=[]
    for i in file_input:
        if ">" in i:
            labels.append(i.replace(">", "").replace('\n', ''))
        else:
            sequences.append(i.replace('\n', ''))
    if label=="labels":
        return labels
    elif label=="sequence":
        return sequences

# Model building
def build_model(file_input, rank):
    # Reads in saved regression model

    unique_sequences=pd.read_csv("unique_sequences_5.csv")
    unique_sequences=unique_sequences.iloc[:,1].values

    input_sequence_list=processor(file_input,label="sequence")
    label_list=processor(file_input, label="labels")
    prediction_output=[]
    significance_output=[]
    taxa_output=[]
    

    for j in input_sequence_list:
        unique_sequences_freq=[]
        unique_sequences_freq.append([j.count(i) for i in unique_sequences])
        pickle_file=open('finalized_model_5_PKJK.pkl','rb')
        model=pickle.load(pickle_file)
        prediction = model.predict(pd.DataFrame(unique_sequences_freq))[0]
        prediction_output.append(prediction)
        significance_output.append(significance(prediction))

        if rank == 'Kingdom':
            taxa_output.append(taxa_extractor(unique_sequences_freq,level=0))
        elif rank == 'Phylum':
            taxa_output.append(taxa_extractor(unique_sequences_freq,level=1))
        elif rank == 'Class':
            taxa_output.append(taxa_extractor(unique_sequences_freq,level=2))
        elif rank == 'Order':
            taxa_output.append(taxa_extractor(unique_sequences_freq,level=3))
        elif rank == 'Family':
            taxa_output.append(taxa_extractor(unique_sequences_freq,level=4))
        elif rank == 'Genus':
            taxa_output.append(taxa_extractor(unique_sequences_freq,level=5))
            
        
    st.header('**Prediction output**')

    prediction_output = pd.concat([ pd.DataFrame(input_sequence_list), pd.DataFrame(prediction_output), pd.DataFrame(taxa_output) , pd.DataFrame(significance_output)  ], axis=1)
    prediction_output.columns = ['Sequence Input','pKJK5','Closest '+ rank,'Greater Than % Baseline Population']
    prediction_output.index=label_list

    st.write(prediction_output)
    st.markdown(filedownload(prediction_output), unsafe_allow_html=True)
 

image = Image.open('logo.png')
st.image(image, use_column_width=False)

# Page title
st.markdown("""
# Plasmid Permissiveness Prediction
This app allows you to predict the plasmid permissiveness in Wastewater Treatement Plants from the partial `16s rRNA` sequence. The prediction is made for three broad-host range plasmids of `pKJK5`, `pB10` and `RP4`.

**Credits**
- App built in `Python` + `Streamlit` by the [Moradigaravand](https://www.birmingham.ac.uk/staff/profiles/cancer-genomic/moradigaravand-danesh.aspx) lab
---
""")

# Sidebar , type=['txt']
with st.sidebar.header('Upload your 16s rRNA data (multifasta format)'):
    uploaded_file = st.sidebar.file_uploader("Upload your input file")
    st.sidebar.markdown("""
[Example input file](https://raw.githubusercontent.com/DaneshMoradigaravand/PlasmidPerm/main/input.fasta)
""")
with st.sidebar.header('2. What level to show'):
    rank = st.selectbox(
    'What level to show?',
    ['Kingdom'	,'Phylum'	,'Class'	,'Order'	,'Family'	,'Genus'])


if st.sidebar.button('Predict'):
    if uploaded_file is not None:
        stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
        uploaded_file_input=stringio.read().splitlines()

        build_model(uploaded_file_input,rank)
    else:
        st.info('Upload input data in the sidebar to start!')
else:
    st.info('Upload input data in the sidebar to start!')