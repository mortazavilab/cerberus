import streamlit as st

def check_multiselect(var, var_name):
    if len(var) > 1:
        st.write(f'More than 1 {var_name} selected. Please choose just 1.')
        var = None
    elif len(var) == 0:
        var = None
    else:
        var = var[0]
    return var
