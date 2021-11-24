import streamlit as st
import pandas as pd
from PIL import Image
import subprocess
import os
import base64
import pickle

from padelpy import padeldescriptor
import glob

# Molecular descriptor calculator


def desc_calc(data):
    # Performs the descriptor calculation

    selection = ["canonical_smiles", "molecule_chembl_id"]
    df3_selection = df3[selection]
    df3_selection.to_csv("molecule.smi", sep="\t", index=False, header=False)
    xml_files = glob.glob("*.xml")
    xml_files.sort()
    fplist = [
        "AtomPairs2Dcount",
        "AtomPairs2D",
        "EState",
        "Extended",
        "Fingerprint",
        "GraphOnly",
        "KlekotaRothcount",
        "KlekotaRoth",
        "MACCS",
        "Pubchem",
        "SubstructureCount",
        "Substructure",
        "descriptors",
    ]
    fp = dict(zip(fplist, xml_files))
    fingerprint = "Pubchem"

    fingerprint_output_file = "".join([fingerprint, ".csv"])  # Pubchem.csv
    fingerprint_descriptortypes = fp[fingerprint]

    padeldescriptor(
        mol_dir="molecule.smi",
        d_file=fingerprint_output_file,
        # descriptortypes='PubchemFingerprint.xml',
        descriptortypes=fingerprint_descriptortypes,
        detectaromaticity=True,
        standardizenitro=True,
        standardizetautomers=True,
        threads=2,
        removesalt=True,
        log=True,
        fingerprints=True,
    )

    os.remove("molecule.smi")


# File download


def filedownload(df):
    csv = df.to_csv(index=False)
    # strings <-> bytes conversions
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href


# Model building


def build_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open("acetylcholinesterase_model.pkl", "rb"))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header("**Prediction output**")
    prediction_output = pd.Series(prediction, name="pIC50")
    molecule_name = pd.Series(load_data[1], name="molecule_name")
    df = pd.concat([molecule_name, prediction_output], axis=1)
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)


# Logo image
image = Image.open("logo.png")

st.image(image, use_column_width=True)

# Page title
st.markdown(
    """
# Bioactivity Prediction App (Acetylcholinesterase)

This app allows you to predict the bioactivity towards inhibting the `Acetylcholinesterase` enzyme. `Acetylcholinesterase` is a drug target for Alzheimer's disease.

**Credits**
- App built in `Python` + `Streamlit`
- Descriptor calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/)
---
"""
)

# Sidebar
with st.sidebar.header("1. Upload your CSV data"):
    uploaded_file = st.sidebar.file_uploader("Upload your input file", type=["txt"])
    st.sidebar.markdown(
        """
[Example input file](https://raw.githubusercontent.com/dataprofessor/bioactivity-prediction-app/main/example_acetylcholinesterase.txt)
"""
    )

if st.sidebar.button("Predict"):
    load_data = pd.read_table(uploaded_file, sep=" ", header=None)
    load_data.to_csv("molecule.smi", sep="\t", header=False, index=False)

    st.header("**Original input data**")
    st.write(load_data)

    with st.spinner("Calculating descriptors..."):
        desc_calc(load_data)

    # Read in calculated descriptors and display the dataframe
    st.header("**Calculated molecular descriptors**")
    desc = pd.read_csv("pubchem.csv")
    st.write(desc)
    st.write(desc.shape)

    # Read descriptor list used in previously built model
    st.header("**Subset of descriptors from previously built models**")
    Xlist = list(pd.read_csv("descriptor_list.csv").columns)
    desc_subset = desc[Xlist]
    st.write(desc_subset)
    st.write(desc_subset.shape)

    # Apply trained model to make prediction on query compounds
    build_model(desc_subset)
else:
    st.info("Upload input data in the sidebar to start!")
