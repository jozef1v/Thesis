from gooey import Gooey, GooeyParser
import pandas as pd
import numpy as np
from cvxopt import matrix, solvers, printing
import xml.etree.ElementTree as ET
import struct
import binascii
from tabulate import tabulate
import warnings
warnings.filterwarnings("error")
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TKAgg')
from pprintpp import pprint
import time
import sys
import os
from scipy.stats import spearmanr, pearsonr, kendalltau, entropy, kstest
from Levenshtein import distance
import re
from colored import stylize, attr, fg
import json
from operator import itemgetter
from mosek import iparam

if len(sys.argv) >= 2:
    if not "--ignore-gooey" in sys.argv:
        sys.argv.append("--ignore-gooey")

if not "--ignore-gooey" in sys.argv and getattr(sys, "frozen", False):
    import win32process, win32gui, win32con
    import os

    pid = os.getpid()

    def get_hwnds_for_pid(pid):
        def callback(hwnd, hwnds):
            if win32gui.IsWindowVisible(hwnd) and win32gui.IsWindowEnabled(hwnd):
                _, found_pid = win32process.GetWindowThreadProcessId(hwnd)
                if found_pid == pid:
                    hwnds.append(hwnd)
            return True

        hwnds = []
        win32gui.EnumWindows(callback, hwnds)
        return hwnds

    open_windows_with_corresponding_pid = get_hwnds_for_pid(pid)
    for hwnd in open_windows_with_corresponding_pid:
        try:
            win32gui.ShowWindow(hwnd, win32con.SW_HIDE)
        except Exception:
            pass


def collect_names(root, place, SS_val, unqs):
    names = []

    # Iterate over variables (CV, MV & DV)
    for count, places in enumerate(place):
        names.append([])

        # Iterate over names of variables (CV, MV & DV)
        for result in root.findall(".//T[@D='" + places + "']/T[@Name]"):
            names[count].append(result.attrib['Name'])

    place[1], place[2] = place[2], place[1]
    varbs = [[n[13:-3] for n in SS_val.Variable.tolist()[:unqs[0]]], \
                [n[13:-3] for n in SS_val.Variable.tolist()[unqs[0]:unqs[0] + unqs[1]]], \
                [n[13:-3] for n in SS_val.Variable.tolist()[unqs[0] + unqs[1]:]]]
    rmvs = ["OP", "PV", "SP"]
    for i, vrb in enumerate(varbs):
        for j, vr in enumerate(vrb):
            if vr[-2:] in rmvs:
                varbs[i][j] = vr[:-2]
    varbs[1], varbs[2] = varbs[2], varbs[1]

    index = []
    for i, valI in enumerate(varbs):
        index.append([])
        for j, valJ in enumerate(valI):
            for k, valK in enumerate(names[i]):
                if valJ in valK:
                    index[i].append(k)
    index[1], index[2] = index[2], index[1]

    return index


def collect_data(root, name, place, properties):
    data = {}

    # Iterate over variables (CV, MV & DV)
    for count, places in enumerate(place):
        data[count] = {}

        # Iterate over names of variables (CV, MV & DV)
        for count2, result in enumerate(root.findall(".//T[@D='" + places + "']/T[@Name]")):
            data[count][count2] = {}
            data[count][count2][name[count]] = result.attrib['Name']

            # Iterate over properties
            for property in properties:
                try:

                    # Iterate over properties of names of variables (CV, MV & DV)
                    for result2 in result.find("T[@Name='" + property + "']"):
                        if str(result2.tag) == 'V':

                            # Convert HEX to FLOAT
                            if result2.text[:2] == '0X':
                                data[count][count2][property] = struct.unpack('!f', bytes.fromhex(result2.text[2:]))[0]

                            # Append FLOAT
                            else:
                                data[count][count2][property] = result2.text

                # Property not found
                except:
                    if (name[count] == 'CV' and property in (properties[0:3] + properties[8:len(properties)])) or \
                        (name[count] == 'MV' and property in (properties[3:6] + properties[8:len(properties)])) or \
                        (name[count] == 'DV' and property in properties[6:10]):
                        data[count][count2][property] = 'No value in URT file'

    return data

# Substitute SS values (XLSX over URT)
def insert_SS(data, CV_excel, MV_excel):
    CV_arr, MV_arr = CV_excel.to_numpy(), MV_excel.to_numpy()
    CV_dic, MV_dic = {}, {}

    # Assign to variables (CV & MV) their steady state values
    for i in range(len(CV_arr)):
        CV_dic[str(CV_arr[i,1] + '_' + CV_arr[i,2])] = CV_arr[i,3]
    for i in range(len(MV_arr)):
        MV_dic[str(MV_arr[i,1] + '_' + MV_arr[i,2])] = MV_arr[i,3]

    # Substitute steady state values (XLSX over URT)
    dics = [CV_dic, MV_dic]
    for i in range(len(dics)):
        for j in range(len(dics[i])):
            for k in dics[i].keys():
                if i == 0 and k == data[i][j]['CV'] or i == 1 and k == data[i][j]['MV']:
                    data[i][j]['SS Value'] = float(dics[i][k])
    return data

# Objective function term weights
def weights(root, data):
    load_data = ['CV scaling value', 'CV low limit error weight', 'CV high limit error weight']
    assign_data = ['Scaling', 'EU_low', 'EU_high']
    assign_data2 = ['Weight_low', 'Weight_high']

    # Iterate over load data
    for i in range(len(load_data)):

        # Iterate over weight data
        for count, result in enumerate(root.findall(".//*[@D='" + load_data[i] + "']")):
            for child in result:

                # Convert HEX to FLOAT
                if "0X" in child.text:
                    data[0][count].update({assign_data[i]: struct.unpack('!f', bytes.fromhex(child.text[2:]))[0]})

                # Append FLOAT
                else:
                    data[0][count].update({assign_data[i]: child.text})

    # Iterate over data
    for i in range(1, len(data)):

        # Calculate & append weights
        for count, result in enumerate(root.findall(".//*[@D='" + load_data[i] + "']")):
            data[0][count].update({assign_data2[i-1]: data[0][count]['Scaling']/data[0][count][assign_data[i]]})

    return data


# Collect gain matrix data
def collect_data_for_gain_matrix(root, name, place, gains):
    data = {}

    # Iterate over CV variables
    for num, CV_name in enumerate(root.findall(".//T[@D='" + place + "']/T[@Name]")):
        data[num] = {}

        # Append CV variable name
        data[num][name] = CV_name.attrib['Name']

        # Iterate over variable (MV & DV) gains
        for gain in gains:
            list = []

            # Iterate over gain values
            for child in CV_name.find('T[@Name="' + gain + '"]/T[@Name="Gain"]'):

                # Gain as N-time value
                try:
                    if type(child.text) == type(None):
                        list.append('N' + str(child.attrib['N']) + 'thisempty')
                    else:
                        list.append('N' + str(child.attrib['N']) + 'this' + str(child.text))

                # Gain as one value
                except:
                    if type(child.text) == type(None):
                        list.append('N1thisempty')
                    else:
                        list.append(str(child.text))

            # Append gains
            data[num][gain] = list

    return data


# Load Model Usage data
def model_usage_vector(root):
    m  = []

    # Model usage data
    for value in root.find('.//T[@Name="Model Usage"]'):

        # Model usage data as N-time value
        try:
            m.extend([value.text]*int(value.attrib['N']))

        # Model usage data as one value
        except:
            m.extend(value.text)

    return m

# Model gain position template
def create_model_usage_array(data_gain_matrix, model_usage_vector_data, data, arg_s):

    # Rewrite data gain matrix to array of proper size
    model_usage_template = create_model_usage_template(data_gain_matrix)

    rows, cols = len(model_usage_template), len(model_usage_template[0])
    
    if arg_s.Time == 'Single':
        MV_size = len(data[1])
    elif arg_s.Time == 'Multiple':
        MV_size = len(data)

    # Iterate over model rows
    for i in range(rows):

        # Iterate over model columns
        for j in range(cols):

            # Non-empty value of type '2'
            if model_usage_template[i][j] != 'empty' and model_usage_vector_data[0] == '2':
                model_usage_template[i][j] = 1.0
                model_usage_vector_data.pop(0)

            # Non-empty value
            elif model_usage_template[i][j] != 'empty' and model_usage_vector_data[0] != '2':
                model_usage_template[i][j] = 0.0
                model_usage_vector_data.pop(0)

            # Empty value
            elif model_usage_template[i][j] == 'empty':  
                model_usage_template[i][j] = 0.0

    # Model usage array
    array = np.array(model_usage_template, dtype= float)
    model_usage_array = array[:, 0:MV_size]

    return model_usage_array

# Rewrite data gain matrix to array of proper size
def create_model_usage_template(data_gain_matrix):
    y = []
    gains = ['MV Gain Delay', 'DV Gain Delay']

    # Iterate over variable (MV & DV) names gains
    for i in data_gain_matrix.keys():
        x = []

        # Iterate over variables (MV & DV) gains
        for gain in gains:

            # Iterate over gain data
            for var in data_gain_matrix[i][gain]:

                # Gain as N-time value
                if var[0] == 'N':
                    splited = var.split('this')
                    x += [splited[1]]*int(splited[0][1::])

                # Gain as one value
                else:
                    x += [var]

        # Append gain
        y.append(x)

    return y

# Create gain matrix
def create_gain_matrix(data_gain_matrix, gains):
    gain_matrix = []

    # Iterate over variables (MV & DV)
    for i in range(2):
        gain_data = []

        # Iterate over variable gain rows
        for j in range(len(data_gain_matrix)):
            row_data = []

            # Iterate over variable gain row data
            for value in data_gain_matrix[j][gains[i]]:

                # Value is N-fold
                if value[0] == 'N':

                    # Original value was empty
                    if value[-1] == 'y':
                        row_data = row_data + [0.0]*int(value.split('this')[0][1:])

                    # Original value was non-empty
                    else:
                        row_data = row_data + [float(value.split('this')[1])]*int(value.split('this')[0][1:])

                # Value is 1-fold
                else:
                    row_data = row_data + [float(value)]

            # Append row data
            gain_data.append(row_data)

        # Append variable gain data
        gain_matrix.append(np.array(gain_data))

    return gain_matrix

# Filter data by status
def filter_data_by_status(data):

    var = {"CV": {}, "MV": {}, "DV": {}}
    var_index = {"CV": [], "MV": [], "DV": []}

    # Iterate over variables (CV, MV & DV)
    for pos, key in enumerate(var.keys()):
        index = 0

        # Iterate over variable data
        for key2, i in enumerate(data[pos]):

            # Specify CV & MV data
            if (key == "CV" and (int(data[pos][key2]['Status']) == 1 or int(data[pos][key2]['Status']) == 2)) or \
                (key == "MV" and (int(data[pos][key2]['Status']) == 1 or int(data[pos][key2]['Status']) == 2 or int(data[pos][key2]['Status']) == 3)):
                data[pos][key2]['index'] = index
                var[key][key2] = data[pos][key2]
                index += 1

            # Specify DV data
            elif key == "DV" and int(data[pos][key2]['Status']) == 1:
                var[key][key2] = data[pos][key2]

            # Append index
            else:
                var_index[key].append(key2)

    return var['CV'], var_index['CV'], var['MV'], var_index['MV'], var['DV'], var_index['DV'], data

# Fill in the missing bounds
def create_missing_bounds(data_CV, data_MV, differ):

    # Set CV bounds
    for i in data_CV.keys():

        # Real high limits
        try:
            data_CV[i]['High Limit'] = struct.unpack('!f', bytes.fromhex(data_CV[i]['High Limit'][2:]))[0]
        except:
            try:
                data_CV[i]['High Limit'] = float(data_CV[i]['High Limit'])
            except:
                try:
                    data_CV[i]['High Limit'] = struct.unpack('!f', bytes.fromhex(data_CV[i]['Ideal High'][2:]))[0]
                except:
                    try:
                        data_CV[i]['High Limit'] = float(data_CV[i]['Ideal High'])
                    except:
                        try:
                            data_CV[i]['High Limit'] = float(data_CV[i]['SS Value']) + abs(float(data_CV[i]['SS Value']) - float(data_CV[i]['Low Limit']))
                        except:
                            data_CV[i]['High Limit'] = float(data_CV[i]['SS Value']) + float(abs(data_CV[i]['SS Value']))*differ
        # data_CV[i]['High Limit'] = 1e-4

        # Real low limits
        try:
            data_CV[i]['Low Limit'] = struct.unpack('!f', bytes.fromhex(data_CV[i]['Low Limit'][2:]))[0]
        except:
            try:
                data_CV[i]['Low Limit'] = float(data_CV[i]['Low Limit'])
            except:
                try:
                    data_CV[i]['Low Limit'] = struct.unpack('!f', bytes.fromhex(data_CV[i]['Ideal Low'][2:]))[0]
                except:
                    try:
                        data_CV[i]['Low Limit'] = float(data_CV[i]['Ideal Low'])
                    except:
                        try:
                            data_CV[i]['Low Limit'] = float(data_CV[i]['SS Value']) - abs(float(data_CV[i]['SS Value']) - float(data_CV[i]['Low Limit']))
                        except:
                            data_CV[i]['Low Limit'] = float(data_CV[i]['SS Value']) - float(abs(data_CV[i]['SS Value']))*differ

        # Ideal high limits
        try:
            data_CV[i]['Ideal High'] = struct.unpack('!f', bytes.fromhex(data_CV[i]['Ideal High'][2:]))[0]
        except:
            try:
                data_CV[i]['Ideal High'] = float(data_CV[i]['Ideal High'])
            except:
                try:
                    data_CV[i]['Ideal High'] = struct.unpack('!f', bytes.fromhex(data_CV[i]['High Limit'][2:]))[0]
                except:
                    try:
                        data_CV[i]['Ideal High'] = float(data_CV[i]['High Limit'])
                    except:
                        try:
                            data_CV[i]['Ideal High'] = float(data_CV[i]['SS Value']) + abs(float(data_CV[i]['SS Value']) - float(data_CV[i]['Low Limit']))
                        except:
                            data_CV[i]['Ideal High'] = float(data_CV[i]['SS Value']) + float(abs(data_CV[i]['SS Value']))*differ

        # Ideal low limits
        try:
            data_CV[i]['Ideal Low'] = struct.unpack('!f', bytes.fromhex(data_CV[i]['Ideal Low'][2:]))[0]
        except:
            try:
                data_CV[i]['Ideal Low'] = float(data_CV[i]['Ideal Low'])
            except:
                try:
                    data_CV[i]['Ideal Low'] = struct.unpack('!f', bytes.fromhex(data_CV[i]['Low Limit'][2:]))[0]
                except:
                    try:
                        data_CV[i]['Ideal Low'] = float(data_CV[i]['Low Limit'])
                    except:
                        try:
                            data_CV[i]['Ideal Low'] = float(data_CV[i]['SS Value']) - abs(float(data_CV[i]['SS Value']) - float(data_CV[i]['Low Limit']))
                        except:
                            data_CV[i]['Ideal Low'] = float(data_CV[i]['SS Value']) - float(abs(data_CV[i]['SS Value']))*differ


    # Set MV bounds
    for i in data_MV.keys():

        # Real high limits
        try:
            data_MV[i]['High Limit'] = struct.unpack('!f', bytes.fromhex(data_MV[i]['High Limit'][2:]))[0]
        except:
            try:
                data_MV[i]['High Limit'] = float(data_MV[i]['High Limit'])
            except:
                try:
                    data_MV[i]['High Limit'] = struct.unpack('!f', bytes.fromhex(data_MV[i]['Ideal High'][2:]))[0]
                except:
                    try:
                        data_MV[i]['High Limit'] = float(data_MV[i]['Ideal High'])
                    except:
                        try:
                            data_MV[i]['High Limit'] = float(data_MV[i]['SS Value']) + abs(float(data_MV[i]['SS Value']) - float(data_MV[i]['Low Limit']))
                        except:
                            data_MV[i]['High Limit'] = float(data_MV[i]['SS Value']) + float(abs(data_MV[i]['SS Value']))*differ

        # Real low limits
        try:
            data_MV[i]['Low Limit'] = struct.unpack('!f', bytes.fromhex(data_MV[i]['Low Limit'][2:]))[0]
        except:
            try:
                data_MV[i]['Low Limit'] = float(data_MV[i]['Low Limit'])
            except:
                try:
                    data_MV[i]['Low Limit'] = struct.unpack('!f', bytes.fromhex(data_MV[i]['Ideal Low'][2:]))[0]
                except:
                    try:
                        data_MV[i]['Low Limit'] = float(data_MV[i]['Ideal Low'])
                    except:
                        try:
                            data_MV[i]['Low Limit'] = float(data_MV[i]['SS Value']) - abs(float(data_MV[i]['SS Value']) - float(data_MV[i]['Low Limit']))
                        except:
                            data_MV[i]['Low Limit'] = float(data_MV[i]['SS Value']) - float(abs(data_MV[i]['SS Value']))*differ

        # Ideal high limits
        try:
            data_MV[i]['Ideal High'] = struct.unpack('!f', bytes.fromhex(data_MV[i]['Ideal High'][2:]))[0]
        except:
            try:
                data_MV[i]['Ideal High'] = float(data_MV[i]['Ideal High'])
            except:
                try:
                    data_MV[i]['Ideal High'] = struct.unpack('!f', bytes.fromhex(data_MV[i]['High Limit'][2:]))[0]
                except:
                    try:
                        data_MV[i]['Ideal High'] = float(data_MV[i]['High Limit'])
                    except:
                        try:
                            data_MV[i]['Ideal High'] = float(data_MV[i]['SS Value']) + abs(float(data_MV[i]['SS Value']) - float(data_MV[i]['Low Limit']))
                        except:
                            data_MV[i]['Ideal High'] = float(data_MV[i]['SS Value']) + float(abs(data_MV[i]['SS Value']))*differ

        # Ideal low limits
        try:
            data_MV[i]['Ideal Low'] = struct.unpack('!f', bytes.fromhex(data_MV[i]['Ideal Low'][2:]))[0]
        except:
            try:
                data_MV[i]['Ideal Low'] = float(data_MV[i]['Ideal Low'])
            except:
                try:
                    data_MV[i]['Ideal Low'] = struct.unpack('!f', bytes.fromhex(data_MV[i]['Low Limit'][2:]))[0]
                except:
                    try:
                        data_MV[i]['Ideal Low'] = float(data_MV[i]['Low Limit'])
                    except:
                        try:
                            data_MV[i]['Ideal Low'] = float(data_MV[i]['SS Value']) - abs(float(data_MV[i]['SS Value']) - float(data_MV[i]['Low Limit']))
                        except:
                            data_MV[i]['Ideal Low'] = float(data_MV[i]['SS Value']) - float(abs(data_MV[i]['SS Value']))*differ

    return data_CV, data_MV

# Adjust the gain matrix to the final form
def filter_gain_matrix(index_CV, index_MV, index_DV, gain_matrix_model, model_usage_array):

    # Modify gain matrix
    gain_matrix_model[0] = gain_matrix_model[0]*model_usage_array

    # Delete unused gain matrix rows
    gain_MV_r = np.delete(gain_matrix_model[0], index_CV, 0)
    gain_DV_r = np.delete(gain_matrix_model[1], index_CV, 0)

    # Delete unused gain matrix columns
    gain_MV = np.delete(gain_MV_r, index_MV, 1)
    gain_DV = np.delete(gain_DV_r, index_DV, 1)

    return gain_MV, gain_DV

# Create P matrix
def create_matrix_P(data_CV, data_MV):

    # P - bez softu, J. Puk
    # P2 - 2-norma^2
    # P4 - 1-norma

    PC, PC2, PC4 = [], [], []

    # CV quadratic penalisation
    for i in data_CV.keys():
        PC.append(float(data_CV[i]['Quadratic Coeff'])**2)
        PC2.append(float(data_CV[i]['Quadratic Coeff'])**2)
        PC4.append(float(data_CV[i]['Quadratic Coeff'])**2)

    # MV quadratic penalisation
    for i in data_MV.keys():
        PC.append(float(data_MV[i]['Quadratic Coeff'])**2)
        PC2.append(float(data_MV[i]['Quadratic Coeff'])**2)
        PC4.append(float(data_MV[i]['Quadratic Coeff'])**2)

    # CV 2-norm^2 & 1-norm quadratic penalization
    for i in data_CV.keys():
        PC2 = PC2 + 2*[(float(data_CV[i]['Quadratic Coeff'])**2 + 0.001)*1e50]
        PC4 = PC4 + 4*[0.0]

    # MV 2-norm^2 & 1-norm quadratic penalization
    for i in data_MV.keys():
        PC2 = PC2 + 2*[(float(data_MV[i]['Quadratic Coeff'])**2 + 0.001)*1e50]
        PC4 = PC4 + 4*[0.0]

    P = matrix(np.diag(PC))
    P2 = matrix(np.diag(PC2))
    P4 = matrix(np.diag(PC4))

    return P, P2, P4

# Create q matrix
def create_matrix_q(data_CV, data_MV):

    try:
        # q - bez softu, J. Puk
        # q2 - 2-norma^2
        # q4 - 1-norma

        q, q2, q4 = [], [], []

        # CV linear penalisation
        for i in data_CV.keys():
            q.append(float(data_CV[i]['Linear Coeff']) - 2.0*float(data_CV[i]['Quadratic Coeff'])**2*float(data_CV[i]['Desired Value']))
            q2.append(float(data_CV[i]['Linear Coeff']) - 2.0*float(data_CV[i]['Quadratic Coeff'])**2*float(data_CV[i]['Desired Value']))
            q4.append(float(data_CV[i]['Linear Coeff']) - 2.0*float(data_CV[i]['Quadratic Coeff'])**2*float(data_CV[i]['Desired Value']))

        # MV linear penalisation
        for i in data_MV.keys():
            q.append(float(data_MV[i]['Linear Coeff']) - 2.0*float(data_MV[i]['Quadratic Coeff'])**2*float(data_MV[i]['Desired Value']))
            q2.append(float(data_MV[i]['Linear Coeff']) - 2.0*float(data_MV[i]['Quadratic Coeff'])**2*float(data_MV[i]['Desired Value']))
            q4.append(float(data_MV[i]['Linear Coeff']) - 2.0*float(data_MV[i]['Quadratic Coeff'])**2*float(data_MV[i]['Desired Value']))

        # CV 2-norm^2 linear penalisation
        for i in data_CV.keys():
            q2 = q2 + 2*[0.0]
            q4 = q4 + 2*[0.0]

        # MV 2-norm^2 linear penalisation
        for i in data_MV.keys():
            q2 = q2 + 2*[0.0]
            q4 = q4 + 2*[0.0]

        # CV 1-norm linear penalisation
        for i in data_CV.keys():
            # ROVNAKÉ AKO PREDOŠLÉ
            q4.append(10**(len(str(max(([float(data_CV[i]['Linear Coeff']) - 2.0*float(data_CV[i]['Quadratic Coeff'])**2*float(data_CV[i]['Desired Value']) for i in data_CV.keys()])))))*1e10)
            q4.append(10**(len(str(max(([float(data_CV[i]['Linear Coeff']) - 2.0*float(data_CV[i]['Quadratic Coeff'])**2*float(data_CV[i]['Desired Value']) for i in data_CV.keys()])))))*1e10)

        # MV 1-norm linear penalisation
        for i in data_MV.keys():            
            # ROVNAKÉ AKO PREDOŠLÉ
            q4.append(10**(len(str(max(([float(data_MV[i]['Linear Coeff']) - 2.0*float(data_MV[i]['Quadratic Coeff'])**2*float(data_MV[i]['Desired Value']) for i in data_MV.keys()])))))*1e10)
            q4.append(10**(len(str(max(([float(data_MV[i]['Linear Coeff']) - 2.0*float(data_MV[i]['Quadratic Coeff'])**2*float(data_MV[i]['Desired Value']) for i in data_MV.keys()])))))*1e10)

        q = matrix(q)
        q2 = matrix(q2)
        q4 = matrix(q4)

    except Exception as e:
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)

    return q, q2, q4

# Create G matrix
def create_matrix_G(data_CV, data_MV):

    # G - bez softu, J. Puk
    # G2 - 2-norma^2
    # G4 - 1-norma

    size_CV, size_MV = len(data_CV), len(data_MV)

    CV_eye, MV_eye = matrix(np.eye(size_CV, dtype = float)), matrix(np.eye(size_MV, dtype = float))
    CV_zeros, MV_zeros = matrix(np.zeros((size_MV, size_CV), dtype = float)), matrix(np.zeros((size_CV, size_MV), dtype = float))
    CV_zeros2, MV_zeros2 = matrix(np.zeros((size_CV, size_CV), dtype = float)), matrix(np.zeros((size_MV, size_MV), dtype = float))

    G = matrix([[CV_eye, -CV_eye, CV_zeros, CV_zeros],
                [MV_zeros, MV_zeros, MV_eye, -MV_eye]])

    G2 = matrix([[CV_eye, -CV_eye, CV_zeros, CV_zeros, CV_zeros2, CV_zeros2, CV_zeros, CV_zeros, CV_zeros, CV_zeros],
                [MV_zeros, MV_zeros, MV_eye, -MV_eye, MV_zeros, MV_zeros, MV_zeros2, MV_zeros2, MV_zeros2, MV_zeros2],
                [-CV_eye, CV_zeros2, CV_zeros, CV_zeros, -CV_eye, CV_zeros2, CV_zeros, CV_zeros, CV_zeros, CV_zeros],
                [CV_zeros2, -CV_eye, CV_zeros, CV_zeros, CV_zeros2, -CV_eye, CV_zeros, CV_zeros, CV_zeros, CV_zeros],
                [MV_zeros, MV_zeros, -MV_eye, MV_zeros2, MV_zeros, MV_zeros, MV_eye, -MV_eye, MV_zeros2, MV_zeros2],
                [MV_zeros, MV_zeros, MV_zeros2, -MV_eye, MV_zeros, MV_zeros, MV_zeros2, MV_zeros2, MV_eye, -MV_eye]])

    G4 = matrix([[CV_eye, -CV_eye, CV_zeros, CV_zeros, CV_zeros2, CV_zeros2, CV_zeros, CV_zeros, CV_zeros, CV_zeros, CV_zeros2, CV_zeros2, CV_zeros2, CV_zeros2, CV_zeros, CV_zeros, CV_zeros, CV_zeros],
                [MV_zeros, MV_zeros, MV_eye, -MV_eye, MV_zeros, MV_zeros, MV_zeros2, MV_zeros2, MV_zeros2, MV_zeros2, MV_zeros, MV_zeros, MV_zeros, MV_zeros, MV_zeros2, MV_zeros2, MV_zeros2, MV_zeros2],
                [-CV_eye, CV_zeros2, CV_zeros, CV_zeros, -CV_eye, CV_zeros2, CV_zeros, CV_zeros, CV_zeros, CV_zeros, CV_eye, -CV_eye, CV_zeros2, CV_zeros2, CV_zeros, CV_zeros, CV_zeros, CV_zeros],
                [CV_zeros2, -CV_eye, CV_zeros, CV_zeros, CV_zeros2, -CV_eye, CV_zeros, CV_zeros, CV_zeros, CV_zeros, CV_zeros2, CV_zeros2, CV_eye, -CV_eye, CV_zeros, CV_zeros, CV_zeros, CV_zeros],
                [MV_zeros, MV_zeros, -MV_eye, MV_zeros2, MV_zeros, MV_zeros, MV_eye, -MV_eye, MV_zeros2, MV_zeros2, MV_zeros, MV_zeros, MV_zeros, MV_zeros, MV_eye, -MV_eye, MV_zeros2, MV_zeros2],
                [MV_zeros, MV_zeros, MV_zeros2, -MV_eye, MV_zeros, MV_zeros, MV_zeros2, MV_zeros2, MV_eye, -MV_eye, MV_zeros, MV_zeros, MV_zeros, MV_zeros, MV_zeros2, MV_zeros2, MV_eye, -MV_eye],
                [CV_zeros2, CV_zeros2, CV_zeros, CV_zeros, CV_zeros2, CV_zeros2, CV_zeros, CV_zeros, CV_zeros, CV_zeros, -CV_eye, -CV_eye, CV_zeros2, CV_zeros2, CV_zeros, CV_zeros, CV_zeros, CV_zeros],
                [CV_zeros2, CV_zeros2, CV_zeros, CV_zeros, CV_zeros2, CV_zeros2, CV_zeros, CV_zeros, CV_zeros, CV_zeros, CV_zeros2, CV_zeros2, -CV_eye, -CV_eye, CV_zeros, CV_zeros, CV_zeros, CV_zeros],
                [MV_zeros, MV_zeros, MV_zeros2, MV_zeros2, MV_zeros, MV_zeros, MV_zeros2, MV_zeros2, MV_zeros2, MV_zeros2, MV_zeros, MV_zeros, MV_zeros, MV_zeros, -MV_eye, -MV_eye, MV_zeros2, MV_zeros2],
                [MV_zeros, MV_zeros, MV_zeros2, MV_zeros2, MV_zeros, MV_zeros, MV_zeros2, MV_zeros2, MV_zeros2, MV_zeros2, MV_zeros, MV_zeros, MV_zeros, MV_zeros, MV_zeros2, MV_zeros2, -MV_eye, -MV_eye]])

    return G, G2, G4

# Create h matrix
def create_matrix_h(data_CV, data_MV, ideal = 0):

    # h - J. Puk
    # h2 - 2-norma^2
    # h3 - bez softu
    # h4 - 1-norma

    limit = [['High Limit', 'Low Limit'], ['Ideal High', 'Ideal Low']]
    h, h2, h3, h4 = [], [], [], []
    data_h = [data_CV, data_MV]

    # CV & MV high & low limits
    for i in range(len(data_h)):

        # CV & MV high limits
        for j in data_h[i].keys():

            try:
                h.append(float(data_h[i][j][limit[ideal][0]]) - float(data_h[i][j]['Delta Soft High Limit']))
                h2.append(float(data_h[i][j][limit[ideal][0]]) - float(data_h[i][j]['Delta Soft High Limit']))
                h4.append(float(data_h[i][j][limit[ideal][0]]) - float(data_h[i][j]['Delta Soft High Limit']))
                h3.append(float(data_h[i][j][limit[ideal][0]]))

            except:
                if data_h[i][j]['SS Value'] > 0:
                    data_h[i][j][limit[ideal][0]] = float(data_h[i][j]['SS Value'])*10
                    h.append(float(data_h[i][j]['SS Value'])*10)
                    h2.append(float(data_h[i][j]['SS Value'])*10)
                    h4.append(float(data_h[i][j]['SS Value'])*10)
                    h3.append(float(data_h[i][j]['SS Value'])*10)

                elif data_h[i][j]['SS Value'] < 0:
                    data_h[i][j][limit[ideal][0]] = -float(data_h[i][j]['SS Value'])*10
                    h.append(-float(data_h[i][j]['SS Value'])*10)
                    h2.append(-float(data_h[i][j]['SS Value'])*10)
                    h4.append(-float(data_h[i][j]['SS Value'])*10)
                    h3.append(-float(data_h[i][j]['SS Value'])*10)

                else:
                    data_h[i][j][limit[ideal][0]] = 10.0
                    h.append(10.0)
                    h2.append(10.0)
                    h4.append(10.0)
                    h3.append(10.0)

        # CV & MV low limits
        for j in data_h[i].keys():

            try:
                h.append(-float(data_h[i][j][limit[ideal][1]]) - float(data_h[i][j]['Delta Soft Low Limit']))
                h2.append(-float(data_h[i][j][limit[ideal][1]]) - float(data_h[i][j]['Delta Soft Low Limit']))
                h4.append(-float(data_h[i][j][limit[ideal][1]]) - float(data_h[i][j]['Delta Soft Low Limit']))
                h3.append(-(float(data_h[i][j][limit[ideal][1]])))

            except:
                if data_h[i][j]['SS Value'] > 0:
                    data_h[i][j][limit[ideal][1]] = -float(data_h[i][j]['SS Value'])*10
                    h.append(float(data_h[i][j]['SS Value'])*10)
                    h2.append(float(data_h[i][j]['SS Value'])*10)
                    h4.append(float(data_h[i][j]['SS Value'])*10)
                    h3.append(float(data_h[i][j]['SS Value'])*10)

                elif data_h[i][j]['SS Value'] < 0:
                    data_h[i][j][limit[ideal][1]] = float(data_h[i][j]['SS Value'])*10
                    h.append(-float(data_h[i][j]['SS Value'])*10)
                    h2.append(-float(data_h[i][j]['SS Value'])*10)
                    h4.append(-float(data_h[i][j]['SS Value'])*10)
                    h3.append(-float(data_h[i][j]['SS Value'])*10)

                else:
                    data_h[i][j][limit[ideal][1]] = -10.0
                    h.append(-10.0)
                    h2.append(-10.0)
                    h4.append(-10.0)
                    h3.append(-10.0)

    # 1-norm & 2-norm HIGH & LOW high & low limits
    for i in range(len(data_h)):

        # 1-norm & 2-norm HIGH high limits
        for j in data_h[i].keys():
            try:

                # High limits only of MVs
                if i == 1:
                    h2.append(float(data_h[i][j]['Delta Soft High Limit']))
                    h4.append(float(data_h[i][j]['Delta Soft High Limit']))

            except:
                if data_h[i][j]['SS Value'] > 0:
                    data_h[i][j][limit[ideal][0]] = float(data_h[i][j]['SS Value'])*10

                    # High limits only of MVs
                    if i == 1:
                        h2.append(float(data_h[i][j]['SS Value'])*10)
                        h4.append(float(data_h[i][j]['SS Value'])*10)

                elif data_h[i][j]['SS Value'] < 0:
                    data_h[i][j][limit[ideal][0]] = -float(data_h[i][j]['SS Value'])*10

                    # High limits only of MVs
                    if i == 1:
                        h2.append(-float(data_h[i][j]['SS Value'])*10)
                        h4.append(-float(data_h[i][j]['SS Value'])*10)

                else:
                    data_h[i][j][limit[ideal][0]] = 10.0

                    # High limits only of MVs
                    if i == 1:
                        h2.append(10.0)
                        h4.append(10.0)

        # 1-norm & 2-norm HIGH low limits
        for j in data_h[i].keys():
            h2.append(0.0)
            h4.append(0.0)

        # 1-norm & 2-norm LOW high limits
        for j in data_h[i].keys():
            try:

                # High limits only of MVs
                if i == 1:
                    h2.append(float(data_h[i][j]['Delta Soft Low Limit']))
                    h4.append(float(data_h[i][j]['Delta Soft Low Limit']))

            except:
                if data_h[i][j]['SS Value'] > 0:
                    data_h[i][j][limit[ideal][1]] = -float(data_h[i][j]['SS Value'])*10

                    # High limits only of MVs
                    if i == 1:
                        h2.append(float(data_h[i][j]['SS Value'])*10)
                        h4.append(float(data_h[i][j]['SS Value'])*10)

                elif data_h[i][j]['SS Value'] < 0:
                    data_h[i][j][limit[ideal][1]] = float(data_h[i][j]['SS Value'])*10

                    # High limits only of MVs
                    if i == 1:
                        h2.append(-float(data_h[i][j]['SS Value'])*10)
                        h4.append(-float(data_h[i][j]['SS Value'])*10)

                else:
                    data_h[i][j][limit[ideal][0]] = -10.0

                    # High limits only of MVs
                    if i == 1:
                        h2.append(-10.0)
                        h4.append(-10.0)

        # 1-norm & 2-norm LOW low limits
        for j in data_h[i].keys():
            h2.append(0.0)
            h4.append(0.0)

    # 1-norm HIGH & LOW high & low limits
    for i in range(len(data_h)):
        for j in data_h[i].keys():
            h4 = h4 + 4*[0.0]

    h = matrix(h)
    h2 = matrix(h2)
    h3 = matrix(h3)
    h4 = matrix(h4)

    return h, h2, h3, h4, data_CV, data_MV

# Create r
def create_matrix_r(data_CV, data_MV):
    r = []
    data = [data_CV, data_MV]

    # Iterate over variables (CV, MV & DV)
    for i in range(len(data)):

        # Iterate over variable data
        for j in data[i].keys():
            r.append(float(data[i][j]['Quadratic Coeff'])**2*float(data[i][j]['Desired Value'])**2)

    # Scalar value
    r = matrix(np.sum(r, axis = 0))

    return r

# Create A matrix
def create_matrix_A(data_CV, gain_MV, data_MV):

    # A - bez softu, J. Puk
    # A2 - 2-norma^2
    # A4 - 1-norma

    CV_size, MV_size = len(data_CV), len(data_MV)
    CV_eye = matrix(np.eye(CV_size, dtype = float))
    MV = matrix(-gain_MV)

    A = matrix([[CV_eye], [MV]])
    A2 = matrix([[CV_eye], [MV], [matrix(np.zeros((CV_size, 2*CV_size), dtype = float))], [matrix(np.zeros((CV_size, 2*MV_size), dtype = float))]])
    A4 = matrix([[CV_eye], [MV], [matrix(np.zeros((CV_size, 4*CV_size), dtype = float))], [matrix(np.zeros((CV_size, 4*MV_size), dtype = float))]])

    return A, A2, A4

# Create b matrix
def create_matrix_b(data_CV, data_MV, gain_MV):
    CV_size = len(data_CV)
    b, CV_steady_state, MV_steady_state = [], [], []

    # CV steady states
    for i in data_CV.keys():
        CV_steady_state.append(data_CV[i]['SS Value'])
    CV_SS = np.array(CV_steady_state, dtype = float)

    # MV steady states
    for i in data_MV.keys():
        MV_steady_state.append(data_MV[i]['SS Value'])
    MV_SS = np.array(MV_steady_state, dtype = float)

    for i in range(CV_size):
        b.append((CV_SS[i] - np.sum(gain_MV[i]*MV_SS)))
    b = matrix(np.array(b, dtype = float))

    return b

# Insert steady state data
def insert_solution_into_data(solution, solution2, solution3, solution4, data_CV, data_MV, prob, round_to):
    bound_type = [' hard', '', ' 2-norm', ' 1-norm']
    data_type = [data_CV, data_MV]
    sol_type = [solution3, solution, solution2, solution4]

    # Iterate over boundary types
    for l, k in enumerate(bound_type):

        # Iterate over variables (CV & MV)
        for m, j in enumerate(data_type):

            # Iterate over variable data
            for h, i in enumerate(j.keys()):
                if prob[0][l] == 0:
                    j[i]['SS Real' + k] = np.around(sol_type[l]['real']['x'][h + m*len(data_CV)], round_to)
                else:
                    j[i]['SS Real' + k] = 0.0
                if prob[1][l] == 0:
                    j[i]['SS Ideal' + k] = np.around(sol_type[l]['ideal']['x'][h + m*len(data_CV)], round_to)
                else:
                    j[i]['SS Ideal' + k] = 0.0

    return data_CV, data_MV

def var_names():
    variables = ['CV',
                 'MV',
                 'DV']

    place = ['Array of CV info structures',
             'Array of MV info structures',
             'Array of DV info structures']

    gains = ['MV Gain Delay',
             'DV Gain Delay']

    properties = ['CV Name',
                  'CV Description',
                  'CV Index',
                  'MV Name',
                  'MV Description',
                  'MV Index',
                  'DV Name',
                  'DV Description',
                  'Status',
                  'Engineering Units',
                  'Low Limit',
                  'High Limit',
                  'Linear Coeff',
                  'Quadratic Coeff',
                  'Desired Value',
                  'Delta Soft High Limit',
                  'Delta Soft Low Limit',
                  'Ideal Low',
                  'Ideal High']

    options_model = ['Real',
                     'Ideal',
                     'Analysis']

    return variables, place, gains, properties, options_model


def swap_data(dataS, position):
    data_new = {}
    for i, I in enumerate(dataS.values()):
        d_type = {}
        for j, J in enumerate(I.values()):
            d_type[position[i][j]] = J
            if i == 0:
                d_type[position[i][j]]['CV Index'] = str(position[i][j]+1)
            elif i == 2:
                d_type[position[i][j]]['MV Index'] = str(position[i][j]+1)
            elif i == 1:
                d_type[position[i][j]]['DV Index'] = str(position[i][j]+1)
        data_new[i] = d_type
    data_new[1], data_new[2] = data_new[2], data_new[1]

    return data_new


def collect_data2(cnf_val, SS_val, argums, root, place):
    lens = [len([value for value in cnf_val['Variable'].tolist() if "_CV_" in value])/len([value for value in SS_val['Variable'].tolist() if "_CV_" in value]), \
    len([value for value in cnf_val['Variable'].tolist() if "_DV_" in value])/len([value for value in SS_val['Variable'].tolist() if "_DV_" in value]), \
    len([value for value in cnf_val['Variable'].tolist() if "_MV_" in value])/len([value for value in SS_val['Variable'].tolist() if "_MV_" in value])]

    unqs = [len([value for value in SS_val['Variable'].tolist() if "_CV_" in value]), \
            len([value for value in SS_val['Variable'].tolist() if "_DV_" in value]), \
            len([value for value in SS_val['Variable'].tolist() if "_MV_" in value])]

    vals = [cnf_val[cnf_val['Variable'].str.contains("_CV_")], \
            cnf_val[cnf_val['Variable'].str.contains("_DV_")], \
            cnf_val[cnf_val['Variable'].str.contains("_MV_")]]

    SSs = [SS_val[SS_val['Variable'].str.contains("_CV_")], \
            SS_val[SS_val['Variable'].str.contains("_DV_")], \
            SS_val[SS_val['Variable'].str.contains("_MV_")]]

    indexes = collect_names(root, place, SS_val, unqs)

    vars = ['CV', 'DV', 'MV']

    data = {}
    for clm in range(1,len(cnf_val.columns)):
        day = {}
        for val_N, val in enumerate(vals):
            vrbs = {}
            for vrb_N in range(unqs[val_N]):
                vrb = {}
                if vars[val_N] != 'DV':
                    vrb['Desired Value'] = val.iloc[int(lens[val_N])*vrb_N, clm]
                    vrb['Delta Soft High Limit'] = val.iloc[int(lens[val_N])*vrb_N + 1, clm]
                    vrb['Delta Soft Low Limit'] = val.iloc[int(lens[val_N])*vrb_N + 2, clm]
                    vrb['High Limit'] = val.iloc[int(lens[val_N])*vrb_N + 3, clm]
                    vrb['Ideal High'] = val.iloc[int(lens[val_N])*vrb_N + 4, clm]
                    vrb['Ideal Low'] = val.iloc[int(lens[val_N])*vrb_N + 5, clm]
                    vrb['Low Limit'] = val.iloc[int(lens[val_N])*vrb_N + 6, clm]
                    vrb['Linear Coeff'] = val.iloc[int(lens[val_N])*vrb_N + 7, clm]
                    vrb['Quadratic Coeff'] = val.iloc[int(lens[val_N])*vrb_N + 8, clm]
                    vrb['Status'] = val.iloc[int(lens[val_N])*vrb_N + 9, clm]
                else:
                    vrb['Status'] = val.iloc[vrb_N, clm]
                vrb[vars[val_N] + ' Index'] = vrb_N + 1
                vrb[vars[val_N] + ' Name'] = val.iloc[int(lens[val_N])*vrb_N, 0].split('.')[0][13:]
                vrb[vars[val_N]] = vrb[vars[val_N] + ' Name']
                vrb['SS Value'] = SSs[val_N].iloc[vrb_N, clm]
                vrbs[vrb_N] = check_data_logic(vrb, argums.Range)
            if vars[val_N] == 'CV':
                day[0] = vrbs
            elif vars[val_N] == 'MV':
                day[1] = vrbs
            else:
                day[2] = vrbs
        data[clm] = swap_data(day, indexes) # day

    return data

def check_data_logic(vrbD, differ):

    try:
        if vrbD['High Limit'] == 0.0:
            if vrbD['High Limit'] - vrbD['Delta Soft High Limit'] <= vrbD['SS Value'] and \
            vrbD['High Limit'] - vrbD['Delta Soft High Limit'] <= vrbD['Low Limit'] + vrbD['Delta Soft Low Limit']:
                if vrbD['Status'] < vrbD['Ideal High'] - vrbD['Delta Soft High Limit'] and \
                vrbD['Low Limit'] + vrbD['Delta Soft Low Limit'] < vrbD['Ideal High'] - vrbD['Delta Soft High Limit']:
                    vrbD['High Limit'] = vrbD['Ideal High']
                elif vrbD['Status'] < vrbD['Status'] + abs(vrbD['Status'] - vrbD['Low Limit']) - vrbD['Delta Soft High Limit'] and \
                vrbD['Low Limit'] + vrbD['Delta Soft Low Limit'] < vrbD['Status'] + abs(vrbD['Status'] - vrbD['Low Limit']) - vrbD['Delta Soft High Limit']:
                    vrbD['High Limit'] = vrbD['Status'] + abs(vrbD['Status'] - vrbD['Low Limit'])
                elif vrbD['Low Limit'] + vrbD['Delta Soft Low Limit'] <= vrbD['SS Value']:
                    vrbD['High Limit'] = vrbD['Status'] + abs(vrbD['Status'])*differ
                else:
                    vrbD['High Limit'] = vrbD['Low Limit'] + vrbD['Delta Soft Low Limit'] + abs(vrbD['Low Limit'] + vrbD['Delta Soft Low Limit'])*differ
            elif vrbD['High Limit'] - vrbD['Delta Soft High Limit'] <= vrbD['SS Value']:
                if vrbD['Status'] < vrbD['Ideal High'] - vrbD['Delta Soft High Limit']:
                    vrbD['High Limit'] = vrbD['Ideal High']
                elif vrbD['Status'] < vrbD['Status'] + abs(vrbD['Status'] - vrbD['Low Limit']) - vrbD['Delta Soft High Limit']:
                    vrbD['High Limit'] = vrbD['Status'] + abs(vrbD['Status'] - vrbD['Low Limit'])
                else:
                    vrbD['High Limit'] = vrbD['Status'] + abs(vrbD['Status'])*differ
            elif vrbD['High Limit'] - vrbD['Delta Soft High Limit'] <= vrbD['Low Limit'] + vrbD['Delta Soft Low Limit']:
                if vrbD['Low Limit'] + vrbD['Delta Soft Low Limit'] < vrbD['Ideal High'] - vrbD['Delta Soft High Limit']:
                    vrbD['High Limit'] = vrbD['Ideal High']
                elif vrbD['Low Limit'] + vrbD['Delta Soft Low Limit'] < vrbD['Status'] + abs(vrbD['Status'] - vrbD['Low Limit']) - vrbD['Delta Soft High Limit']:
                    vrbD['High Limit'] = vrbD['Status'] + abs(vrbD['Status'] - vrbD['Low Limit'])
                else:
                    vrbD['High Limit'] = vrbD['Low Limit'] + vrbD['Delta Soft Low Limit'] + abs(vrbD['Low Limit'] + vrbD['Delta Soft Low Limit'])*differ
        if vrbD['Low Limit'] == 0.0:
            if vrbD['SS Value'] <= vrbD['Low Limit'] + vrbD['Delta Soft High Limit'] and \
            vrbD['High Limit'] - vrbD['Delta Soft High Limit'] <= vrbD['Low Limit'] + vrbD['Delta Soft Low Limit']:
                if vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit'] < vrbD['Status'] and \
                vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit'] < vrbD['High Limit'] - vrbD['Delta Soft High Limit']:
                    vrbD['Low Limit'] = vrbD['Ideal Low']
                elif vrbD['Status'] - abs(vrbD['Status'] - vrbD['High Limit']) + vrbD['Delta Soft Low Limit'] < vrbD['Status'] and \
                vrbD['Status'] - abs(vrbD['Status'] - vrbD['High Limit']) + vrbD['Delta Soft Low Limit'] < vrbD['High Limit'] + vrbD['Delta Soft High Limit']:
                    vrbD['Low Limit'] = vrbD['Status'] - abs(vrbD['Status'] - vrbD['High Limit'])
                elif vrbD['SS Value'] <= vrbD['High Limit'] - vrbD['Delta Soft High Limit']:
                    vrbD['Low Limit'] = vrbD['Status'] - abs(vrbD['Status'])*differ
                else:
                    vrbD['Low Limit'] = vrbD['High Limit'] - vrbD['Delta Soft High Limit'] - abs(vrbD['High Limit'] - vrbD['Delta Soft High Limit'])*differ
            elif vrbD['SS Value'] <= vrbD['Low Limit'] + vrbD['Delta Soft Low Limit']:
                if vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit'] < vrbD['Status']:
                    vrbD['Low Limit'] = vrbD['Ideal Low']
                elif vrbD['Status'] - abs(vrbD['Status'] - vrbD['High Limit']) + vrbD['Delta Soft Low Limit'] < vrbD['Status']:
                    vrbD['Low Limit'] = vrbD['Status'] - abs(vrbD['Status'] - vrbD['High Limit'])
                else:
                    vrbD['Low Limit'] = vrbD['Status'] - abs(vrbD['Status'])*differ
            elif vrbD['High Limit'] - vrbD['Delta Soft High Limit'] <= vrbD['Low Limit'] + vrbD['Delta Soft Low Limit']:
                if vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit'] < vrbD['High Limit'] - vrbD['Delta Soft High Limit']:
                    vrbD['Low Limit'] = vrbD['Ideal Low']
                elif vrbD['Status'] - abs(vrbD['Status'] - vrbD['High Limit']) + vrbD['Delta Soft Low Limit'] < vrbD['High Limit'] + vrbD['Delta Soft High Limit']:
                    vrbD['Low Limit'] = vrbD['Status'] - abs(vrbD['Status'] - vrbD['High Limit'])
                else:
                    vrbD['Low Limit'] = vrbD['High Limit'] - vrbD['Delta Soft High Limit'] - abs(vrbD['High Limit'] - vrbD['Delta Soft High Limit'])*differ
        if vrbD['Ideal High'] == 0.0:
            if vrbD['Ideal High'] - vrbD['Delta Soft High Limit'] <= vrbD['SS Value'] and \
            vrbD['Ideal High'] - vrbD['Delta Soft High Limit'] <= vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit']:
                if vrbD['Status'] < vrbD['High Limit'] - vrbD['Delta Soft High Limit'] and \
                vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit'] < vrbD['High Limit'] - vrbD['Delta Soft High Limit']:
                    vrbD['Ideal High'] = vrbD['High Limit']
                elif vrbD['Status'] < vrbD['Status'] + abs(vrbD['Status'] - vrbD['Ideal Low']) - vrbD['Delta Soft High Limit'] and \
                vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit'] < vrbD['Status'] + abs(vrbD['Status'] - vrbD['Ideal Low']) - vrbD['Delta Soft High Limit']:
                    vrbD['Ideal High'] = vrbD['Status'] + abs(vrbD['Status'] - vrbD['Ideal Low'])
                elif vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit'] <= vrbD['SS Value']:
                    vrbD['Ideal High'] = vrbD['Status'] + abs(vrbD['Status'])*differ
                else:
                    vrbD['Ideal High'] = vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit'] + abs(vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit'])*differ
            elif vrbD['Ideal High'] - vrbD['Delta Soft High Limit'] <= vrbD['SS Value']:
                if vrbD['Status'] < vrbD['High Limit'] - vrbD['Delta Soft High Limit']:
                    vrbD['Ideal High'] = vrbD['High Limit']
                elif vrbD['Status'] < vrbD['Status'] + abs(vrbD['Status'] - vrbD['Ideal Low']) - vrbD['Delta Soft High Limit']:
                    vrbD['Ideal High'] = vrbD['Status'] + abs(vrbD['Status'] - vrbD['Ideal Low'])
                else:
                    vrbD['Ideal High'] = vrbD['Status'] + abs(vrbD['Status'])*differ
            elif vrbD['Ideal High'] - vrbD['Delta Soft High Limit'] <= vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit']:
                if vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit'] < vrbD['High Limit'] - vrbD['Delta Soft High Limit']:
                    vrbD['Ideal High'] = vrbD['High Limit']
                elif vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit'] < vrbD['Status'] + abs(vrbD['Status'] - vrbD['Ideal Low']) - vrbD['Delta Soft High Limit']:
                    vrbD['Ideal High'] = vrbD['Status'] + abs(vrbD['Status'] - vrbD['Ideal Low'])
                else:
                    vrbD['Ideal High'] = vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit'] + abs(vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit'])*differ
        if vrbD['Ideal Low'] == 0.0:
            if vrbD['SS Value'] <= vrbD['Ideal Low'] + vrbD['Delta Soft High Limit'] and \
            vrbD['Ideal High'] - vrbD['Delta Soft High Limit'] <= vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit']:
                if vrbD['Low Limit'] + vrbD['Delta Soft Low Limit'] < vrbD['Status'] and \
                vrbD['Low Limit'] + vrbD['Delta Soft Low Limit'] < vrbD['Ideal High'] - vrbD['Delta Soft High Limit']:
                    vrbD['Ideal Low'] = vrbD['Low Limit']
                elif vrbD['Status'] - abs(vrbD['Status'] - vrbD['Ideal High']) + vrbD['Delta Soft Low Limit'] < vrbD['Status'] and \
                vrbD['Status'] - abs(vrbD['Status'] - vrbD['Ideal High']) + vrbD['Delta Soft Low Limit'] < vrbD['Ideal High'] + vrbD['Delta Soft High Limit']:
                    vrbD['Ideal Low'] = vrbD['Status'] - abs(vrbD['Status'] - vrbD['Ideal High'])
                elif vrbD['SS Value'] <= vrbD['Ideal High'] - vrbD['Delta Soft High Limit']:
                    vrbD['Ideal Low'] = vrbD['Status'] - abs(vrbD['Status'])*differ
                else:
                    vrbD['Ideal Low'] = vrbD['Ideal High'] - vrbD['Delta Soft High Limit'] - abs(vrbD['Ideal High'] - vrbD['Delta Soft High Limit'])*differ
            elif vrbD['SS Value'] <= vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit']:
                if vrbD['Low Limit'] + vrbD['Delta Soft Low Limit'] < vrbD['Status']:
                    vrbD['Ideal Low'] = vrbD['Low Limit']
                elif vrbD['Status'] - abs(vrbD['Status'] - vrbD['Ideal High']) + vrbD['Delta Soft Low Limit'] < vrbD['Status']:
                    vrbD['Ideal Low'] = vrbD['Status'] - abs(vrbD['Status'] - vrbD['Ideal High'])
                else:
                    vrbD['Ideal Low'] = vrbD['Status'] - abs(vrbD['Status'])*differ
            elif vrbD['Ideal High'] - vrbD['Delta Soft High Limit'] <= vrbD['Ideal Low'] + vrbD['Delta Soft Low Limit']:
                if vrbD['Low Limit'] + vrbD['Delta Soft Low Limit'] < vrbD['Ideal High'] - vrbD['Delta Soft High Limit']:
                    vrbD['Ideal Low'] = vrbD['Low Limit']
                elif vrbD['Status'] - abs(vrbD['Status'] - vrbD['Ideal High']) + vrbD['Delta Soft Low Limit'] < vrbD['Ideal High'] + vrbD['Delta Soft High Limit']:
                    vrbD['Ideal Low'] = vrbD['Status'] - abs(vrbD['Status'] - vrbD['Ideal High'])
                else:
                    vrbD['Ideal Low'] = vrbD['Ideal High'] - vrbD['Delta Soft High Limit'] - abs(vrbD['Ideal High'] - vrbD['Delta Soft High Limit'])*differ
    except:
        pass

    return vrbD

def execute(args_data):

    # Create global variables
    global solution, solution2, solution3, solution4, M, M2, M3, M4, data_CV, data_MV

    try:
        # Load files
        fileExcel = re.search('SteadyStateBFCCU2C1.xlsx', args_data.xlsx).group(0)
        fileURT = re.search('BFCCU2C1_data(.+?).urt', args_data.urt).group(0)

        # Import data from Excel
        CV_excel = pd.read_excel(fileExcel, 'CV')
        MV_excel = pd.read_excel(fileExcel, 'MV')

        # Import data from URT
        tree = ET.parse(fileURT)
        rootET = tree.getroot()

        # Load variables
        variables, place, gains, properties, _ = var_names()

        # Collect CVs, MVs & DVs data
        data = collect_data(rootET, variables, place, properties)

        # Substitute SS values (XLSX over URT)
        data = insert_SS(data, CV_excel, MV_excel)

        # Objective function term weights
        data = weights(rootET, data)

        # Collect gain matrix data
        data_gain_matrix = collect_data_for_gain_matrix(rootET, 'CV', 'Array of CV info structures', gains)

        # Load Model Usage data
        model_usage_vector_data = model_usage_vector(rootET)

        # Model gain position template
        model_usage_array = create_model_usage_array(data_gain_matrix, model_usage_vector_data, data, args_data)

        # Load Model Usage data
        model_usage_vector_data = model_usage_vector(rootET)

        # Create gain matrix
        gain_matrix_model = create_gain_matrix(data_gain_matrix, gains)

        # Filter data by status
        data_CV, index_CV, data_MV, index_MV, data_DV, index_DV, data = filter_data_by_status(data)

        # Fill in the missing bounds
        data_CV, data_MV = create_missing_bounds(data_CV, data_MV, args_data.Range)

        # Adjust the gain matrix to the final form
        gain_MV, gain_DV = filter_gain_matrix(index_CV, index_MV, index_DV, gain_matrix_model, model_usage_array)

        # Column separators
        len_CV, len_MV = len(data_CV), len(data_MV)
        ciarkyCV, ciarkyMV, ciarkyCV2, ciarkyMV2 = ["|"]*len_CV, ["|"]*len_MV, ["||"]*len_CV, ["||"]*len_MV

        # High & low limits
        HCV = [data_CV[i]['High Limit'] for i in data_CV]
        LCV = [data_CV[i]['Low Limit'] for i in data_CV]
        HMV = [data_MV[i]['High Limit'] for i in data_MV]
        LMV = [data_MV[i]['Low Limit'] for i in data_MV]
        iHCV = [data_CV[i]['Ideal High'] for i in data_CV]
        iLCV = [data_CV[i]['Ideal Low'] for i in data_CV]
        iHMV = [data_MV[i]['Ideal High'] for i in data_MV]
        iLMV = [data_MV[i]['Ideal Low'] for i in data_MV]
        HLni = [[HCV, LCV, HMV, LMV], [iHCV, iLCV, iHMV, iLMV]]

        # Delta soft high & low limits
        dHCV = [float(data_CV[i]['Delta Soft High Limit']) for i in data_CV]
        dLCV = [float(data_CV[i]['Delta Soft Low Limit']) for i in data_CV]
        dHMV = [float(data_MV[i]['Delta Soft High Limit']) for i in data_MV]
        dLMV = [float(data_MV[i]['Delta Soft Low Limit']) for i in data_MV]
        dHLni = [[(np.array(HCV) - np.array(dHCV)).tolist(), (np.array(LCV) + np.array(dLCV)).tolist(), (np.array(HMV) - np.array(dHMV)).tolist(), \
                (np.array(LMV) + np.array(dLMV)).tolist()], [(np.array(iHCV) - np.array(dHCV)).tolist(), (np.array(iLCV) + np.array(dLCV)).tolist(), \
                (np.array(iHMV) - np.array(dHMV)).tolist(), (np.array(iLMV) + np.array(dLMV)).tolist()]]

        # Matrices (P, q, r, G, h - ideal & real, A, b)
        P, P2, P4 = create_matrix_P(data_CV, data_MV) 
        q, q2, q4 = create_matrix_q(data_CV, data_MV)
        G, G2, G4 = create_matrix_G(data_CV, data_MV)
        h, h2, h3, h4, data_CV, data_MV = create_matrix_h(data_CV, data_MV)
        h_ideal, h_ideal2, h_ideal3, h_ideal4, data_CV, data_MV = create_matrix_h(data_CV, data_MV, ideal = 1)
        A, A2, A4 = create_matrix_A(data_CV, gain_MV, data_MV)
        b = create_matrix_b(data_CV, data_MV, gain_MV)
        r = create_matrix_r(data_CV, data_MV)

        #   - PUK
        # 2 - 2-norm^2
        # 3 - HARD
        # 4 - 1-norm

        # Matrices dictionaries
        M = {'P': P, 'q': q, 'G': G, 'h': h, 'A': A, 'b': b, 'r': r, 'h_ideal': h_ideal}
        M2 = {'P': P2, 'q': q2, 'G': G2, 'h': h2, 'A': A2, 'b': b, 'r': r, 'h_ideal': h_ideal2}
        M3 = {'P': P, 'q': q, 'G': G, 'h': h3, 'A': A, 'b': b, 'r': r, 'h_ideal': h_ideal3}
        M4 = {'P': P4, 'q': q4, 'G': G4, 'h': h4, 'A': A4, 'b': b, 'r': r, 'h_ideal': h_ideal4}

        # Solve optimization problems (optimal steady states) - real models
        prob_R = [0,0,0,0]
        t_BEZ_SOFTU = time.time();
        try:
            sol3 = solvers.qp(2*P, q, G, h3, A, b, solver = 'mosek')
        except:
            sol3 = 0
            prob_R[0] = 1
        elapsed_BEZ_SOFTU = time.time() - t_BEZ_SOFTU;
        t_PUK_POVODNE = time.time();
        try:
            sol = solvers.qp(2*P, q, G, h, A, b, solver = 'mosek')
        except:
            sol = 0
            prob_R[1] = 1
        elapsed_PUK_POVODNE = time.time() - t_PUK_POVODNE;
        t_2_NORMA = time.time();
        try:
            sol2 = solvers.qp(2*P2, q2, G2, h2, A2, b, solver = 'mosek')
        except:
            sol2 = 0
            prob_R[2] = 1
        elapsed_2_NORMA = time.time() - t_2_NORMA;
        t_1_NORMA = time.time();
        try:
            sol4 = solvers.qp(2*P4, q4, G4, h4, A4, b, solver = 'mosek')
        except:
            sol4 = 0
            prob_R[3] = 1
        elapsed_1_NORMA = time.time() - t_1_NORMA;

        # Solve optimization problems (optimal steady states) - ideal models
        prob_I = [0,0,0,0]
        try:
            sol_ideal3 = solvers.qp(2*P, q, G, h_ideal3, A, b, solver = 'mosek')
        except:
            sol_ideal3 = 0
            prob_I[0] = 1
        try:
            sol_ideal = solvers.qp(2*P, q, G, h_ideal, A, b, solver = 'mosek')
        except:
            sol_ideal = 0
            prob_I[1] = 1
        try:
            sol_ideal2 = solvers.qp(2*P2, q2, G2, h_ideal2, A2, b, solver = 'mosek')
        except:
            sol_ideal2 = 0
            prob_I[2] = 1
        try:
            sol_ideal4 = solvers.qp(2*P4, q4, G4, h_ideal4, A4, b, solver = 'mosek')
        except:
            sol_ideal4 = 0
            prob_I[3] = 1

        print_type = ['\n\n============\n REAL MODEL\n============', '\n\n\n\n=============\n IDEAL MODEL\n=============']
        model_type = ['_real', '_ideal']
        prob = [prob_R, prob_I]
        model_data = [[sol3, sol, sol2, sol4], [sol_ideal3, sol_ideal, sol_ideal2, sol_ideal4]]

        for j, i in enumerate(model_data):
            sol_3, sol_, sol_2, sol_4 = i[0], i[1], i[2], i[3]

            # Model type
            if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                print(print_type[j])

            # Display limit data - command line
            cv_val = {'Variable': ['CV ' + str(data_CV[i]['CV Index']) for i in data_CV], 'Name': [data_CV[i]['CV Name'] for i in data_CV],
                    '||': ciarkyCV2, 'LOW HARD BOUND': HLni[j][1],
                    'LOW SOFT BOUND': dHLni[j][1], '[': ciarkyCV,
                    'HARD': sol_3['x'][:len(data_CV)] if prob[j][0] == 0 else ["--"]*len(data_CV),
                    'TIGHT HARD': sol_['x'][:len(data_CV)] if prob[j][1] == 0 else ["--"]*len(data_CV),
                    '2-norm^2': sol_2['x'][:len(data_CV)] if prob[j][2] == 0 else ["--"]*len(data_CV),
                    '1-norm': sol_4['x'][:len(data_CV)] if prob[j][3] == 0 else ["--"]*len(data_CV),
                    ']': ciarkyCV, 'HIGH SOFT BOUND': dHLni[j][0], 'HIGH HARD BOUND': HLni[j][0]}
            mv_val = {'Variable': ['MV ' + str(data_MV[i]['MV Index']) for i in data_MV], 'Name': [data_MV[i]['MV Name'] for i in data_MV],
                    '||': ciarkyMV2, 'LOW HARD BOUND': HLni[j][3],
                    'LOW SOFT BOUND': dHLni[j][3], '[': ciarkyMV,
                    'HARD': sol_3['x'][len(data_CV):len(data_CV) + len(data_MV)] if prob[j][0] == 0 else ["--"]*len(data_MV),
                    'TIGHT HARD': sol_['x'][len(data_CV):len(data_CV) + len(data_MV)] if prob[j][1] == 0 else ["--"]*len(data_MV),
                    '2-norm^2': sol_2['x'][len(data_CV):len(data_CV) + len(data_MV)] if prob[j][2] == 0 else ["--"]*len(data_MV),
                    '1-norm': sol_4['x'][len(data_CV):len(data_CV) + len(data_MV)] if prob[j][3] == 0 else ["--"]*len(data_MV),
                    ']': ciarkyMV, 'HIGH SOFT BOUND': dHLni[j][2], 'HIGH HARD BOUND': HLni[j][2]}

            # Display limit data - .txt file
            with open("limits_" + fileURT.split('.')[0] + model_type[j] + ".txt", "w") as f:
                f.write('\n==============\nCV CONSTRAINTS\n==============\n\n')
                f.write(tabulate(pd.DataFrame(data = cv_val), headers = 'keys', showindex = False, floatfmt = ".2f"))
                f.write('\n\n\n==============\nMV CONSTRAINTS\n==============\n\n')
                f.write(tabulate(pd.DataFrame(data = mv_val), headers = 'keys', showindex = False, floatfmt = ".2f"))
            if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Solve is True:
                os.system('start "" /max ' + 'limits_' + fileURT.split('.')[0] + model_type[j] + '.txt')
            if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                print()
                print("\nCV optimal\n==========\n")
                print(tabulate(cv_val, headers = 'keys', floatfmt = ".2f"))
                print()
                print("\nMV optimal\n==========\n")
                print(tabulate(mv_val, headers = 'keys', floatfmt = ".2f"))
                print()

            # Display objective function data
            obj = {'HARD': sol_3['x'].trans()*P*sol_3['x'] + q.trans()*sol_3['x'] + r if prob[j][0] == 0 else "----",
                'TIGHT HARD': sol_['x'].trans()*P*sol_['x'] + q.trans()*sol_['x'] + r if prob[j][0] == 0 else "----",
                '2-norm^2 WITH e': sol_2['x'].trans()*P2*sol_2['x'] + q2.trans()*sol_2['x'] + r if prob[j][0] == 0 else "----",
                '2-norm^2 WITHOUT e': sol_2['x'][:len(data_CV) + len(data_MV)].trans()*P2[:len(data_CV) + len(data_MV),:len(data_CV) + len(data_MV)]*sol_2['x'][:len(data_CV) + len(data_MV)] + q2[:len(data_CV) + len(data_MV)].trans()*sol_2['x'][:len(data_CV) + len(data_MV)] + r if prob[j][0] == 0 else "----",
                '1-norm WITH e': sol_4['x'].trans()*P4*sol_4['x'] + q4.trans()*sol_4['x'] + r if prob[j][0] == 0 else "----",
                '1-norm WITHOUT e': sol_4['x'][:len(data_CV) + len(data_MV)].trans()*P4[:len(data_CV) + len(data_MV),:len(data_CV) + len(data_MV)]*sol_4['x'][:len(data_CV) + len(data_MV)] + q4[:len(data_CV) + len(data_MV)].trans()*sol_4['x'][:len(data_CV) + len(data_MV)] + r if prob[j][0] == 0 else "----"}
            if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                print()
                print("\nObjective function value\n========================\n")
                print(tabulate(obj, headers = 'keys', floatfmt = ".4E"))
                print()

            # Display relative objective function data
            rozdiel = {'HARD': (obj['HARD'] - obj['HARD'])/abs(obj['HARD'])*100 if prob[j][0] == 0 else "----",
                    'TIGHT HARD': (obj['TIGHT HARD'] - obj['HARD'])/abs(obj['HARD'])*100 if prob[j][0] == 0 else "----",
                    '2-norm^2 WITH e': (obj['2-norm^2 WITH e'] - obj['HARD'])/abs(obj['HARD'])*100 if prob[j][0] == 0 else "----",
                    '2-norm^2 WITHOUT e': (obj['2-norm^2 WITHOUT e'] - obj['HARD'])/abs(obj['HARD'])*100 if prob[j][0] == 0 else "----",
                    '1-norm WITH e': (obj['1-norm WITH e'] - obj['HARD'])/abs(obj['HARD'])*100 if prob[j][0] == 0 else "----",
                    '1-norm WITHOUT e': (obj['1-norm WITHOUT e'] - obj['HARD'])/abs(obj['HARD'])*100 if prob[j][0] == 0 else "----"}
            if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                print()
                print("\nRelative objective function value (against HARD in %)\n=====================================================\n")
                print(tabulate(rozdiel, headers = 'keys', floatfmt = ".4E"))
                print()

            # Display objective function data - .txt file
            with open("objective_" + fileURT.split('.')[0] + model_type[j] + ".txt", "w") as f:
                f.write('\n===================\nOBJECTIVE FUNCTIONS\n===================\n\n')
                f.write(tabulate(pd.DataFrame(data = obj), headers = 'keys', showindex = False, floatfmt = ".4E"))
                f.write('\n\n\n=====================================================\nRELATIVE OBJECTIVE FUNCTIONS (against BEZ SOFTU IN %)\n=====================================================\n\n')
                f.write(tabulate(pd.DataFrame(data = rozdiel), headers = 'keys', showindex = False, floatfmt = ".4E"))
            if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Solve is True:
                os.system('start "" /max ' + 'objective_' + fileURT.split('.')[0] + model_type[j] + '.txt')

            # Display 2-norm^2 epsilon values
            epsy1 = {'Variable': ['CV ' + str(data_CV[i]['CV Index']) for i in data_CV], 'Name': [data_CV[i]['CV Name'] for i in data_CV],
                    '||': ciarkyCV2, '0': [0]*len_CV, '[': ciarkyCV,
                    '2-norm^2 e_CV,H': sol_2['x'][len(data_CV) + len(data_MV):2*len(data_CV) + len(data_MV)] if prob[j][2] == 0 else ["--"]*len(data_CV),
                    ']': ciarkyCV, 'INF': ['inf']*len_CV}
            epsy2 = {'Variable': ['CV ' + str(data_CV[i]['CV Index']) for i in data_CV], 'Name': [data_CV[i]['CV Name'] for i in data_CV],
                    '||': ciarkyCV2, '0': [0]*len_CV, '[': ciarkyCV,
                    '2-norm^2 e_CV,L': sol_2['x'][2*len(data_CV) + len(data_MV):3*len(data_CV) + len(data_MV)] if prob[j][2] == 0 else ["--"]*len(data_CV),
                    ']': ciarkyCV, 'INF': ['inf']*len_CV}
            epsy3 = {'Variable': ['MV ' + str(data_MV[i]['MV Index']) for i in data_MV], 'Name': [data_MV[i]['MV Name'] for i in data_MV],
                    '||': ciarkyMV2, '0': [0]*len_MV, '[': ciarkyMV,
                    '2-norm^2 e_MV,H': sol_2['x'][3*len(data_CV) + len(data_MV):3*len(data_CV) + 2*len(data_MV)] if prob[j][2] == 0 else ["--"]*len(data_MV),
                    ']': ciarkyMV, 'deltaMV_H': dHMV}
            epsy4 = {'Variable': ['MV ' + str(data_MV[i]['MV Index']) for i in data_MV], 'Name': [data_MV[i]['MV Name'] for i in data_MV],
                    '||': ciarkyMV2, '0': [0]*len_MV, '[': ciarkyMV,
                    '2-norm^2 e_MV,L': sol_2['x'][3*len(data_CV) + 2*len(data_MV):3*len(data_CV) + 3*len(data_MV)] if prob[j][2] == 0 else ["--"]*len(data_MV),
                    ']': ciarkyMV, 'deltaMV_L': dLMV}

            # Display 2-norm^2 epsilon values - .txt file
            with open("epsilons_2_norm_" + fileURT.split('.')[0] + model_type[j] + ".txt", "w") as f:
                f.write('\n=======================\n2-norm^2 epsilon VALUES\n=======================\n\n')
                f.write('\n(0 <= e_CV_H)\n\n')
                f.write(tabulate(pd.DataFrame(data = epsy1), headers = 'keys', showindex = False, floatfmt = ".4E"))
                f.write('\n\n\n(0 <= e_CV_L)\n\n')
                f.write(tabulate(pd.DataFrame(data = epsy2), headers = 'keys', showindex = False, floatfmt = ".4E"))
                f.write('\n\n\n(0 <= e_MV_H <= deltaMV_H)\n\n')
                f.write(tabulate(pd.DataFrame(data = epsy3), headers = 'keys', showindex = False, floatfmt = ".4E"))
                f.write('\n\n\n(0 <= e_MV_L <= deltaMV_L)\n\n')
                f.write(tabulate(pd.DataFrame(data = epsy4), headers = 'keys', showindex = False, floatfmt = ".4E"))
            if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Solve is True:
                os.system('start "" /max ' + 'epsilons_2_norm_' + fileURT.split('.')[0] + model_type[j] + '.txt')
            if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                print(), print("\nEpsilons - 2-norm^2\n===================\n")
                print(), print("(0 <= e_CV_H)\n")
                print(tabulate(epsy1, headers = 'keys', floatfmt = ".4E"))
                print(), print("(0 <= e_CV_L)\n")
                print(tabulate(epsy2, headers = 'keys', floatfmt = ".4E"))
                print(), print("(0 <= e_MV_H <= ΔMV_H)\n")
                print(tabulate(epsy3, headers = 'keys', floatfmt = ".4E"))
                print(), print("(0 <= e_MV_L <= ΔMV_L)\n")
                print(tabulate(epsy4, headers = 'keys', floatfmt = ".4E"))
                print()

            # Display 1-norm epsilon values - .txt file
            epsy1 = {'Variable': ['CV ' + str(data_CV[i]['CV Index']) for i in data_CV], 'Name': [data_CV[i]['CV Name'] for i in data_CV],
                    '||': ciarkyCV2, '-1-norm E_CV,H': sol_4['x'][3*len(data_CV) + 3*len(data_MV):4*len(data_CV) + 3*len(data_MV)] if prob[j][3] == 0 else ["--"]*len(data_CV),
                    '<': ciarkyCV, '0': [0]*len_CV, '[': ciarkyCV,
                    '1-norm e_CV,H': sol_4['x'][len(data_CV) + len(data_MV):2*len(data_CV) + len(data_MV)] if prob[j][3] == 0 else ["--"]*len(data_CV),
                    ']': ciarkyCV, 'INF': ['inf']*len_CV, '>': ciarkyCV,
                    '1-norm E_CV,H': -sol_4['x'][3*len(data_CV) + 3*len(data_MV):4*len(data_CV) + 3*len(data_MV)] if prob[j][3] == 0 else ["--"]*len(data_CV)}
            epsy2 = {'Variable': ['CV ' + str(data_CV[i]['CV Index']) for i in data_CV], 'Name': [data_CV[i]['CV Name'] for i in data_CV],
                    '||': ciarkyCV2, '-1-norm E_CV,L': sol_4['x'][4*len(data_CV) + 3*len(data_MV):5*len(data_CV) + 3*len(data_MV)] if prob[j][3] == 0 else ["--"]*len(data_CV),
                    '<': ciarkyCV, '0': [0]*len_CV, '[': ciarkyCV,
                    '1-norm e_CV,L': sol_4['x'][2*len(data_CV) + len(data_MV):3*len(data_CV) + len(data_MV)] if prob[j][3] == 0 else ["--"]*len(data_CV),
                    ']': ciarkyCV, 'INF': ['inf']*len_CV, '>': ciarkyCV,
                    '1-norm E_CV,L': -sol_4['x'][4*len(data_CV) + 3*len(data_MV):5*len(data_CV) + 3*len(data_MV)] if prob[j][3] == 0 else ["--"]*len(data_CV)}
            epsy3 = {'Variable': ['MV ' + str(data_MV[i]['MV Index']) for i in data_MV], 'Name': [data_MV[i]['MV Name'] for i in data_MV],
                    '||': ciarkyMV2, '-1-norm E_MV,H': -sol_4['x'][5*len(data_CV) + 3*len(data_MV):5*len(data_CV) + 4*len(data_MV)] if prob[j][3] == 0 else ["--"]*len(data_MV),
                    '<': ciarkyMV, '0': [0]*len_MV, '[': ciarkyMV,
                    '1-norm e_MV,H': sol_4['x'][3*len(data_CV) + len(data_MV):3*len(data_CV) + 2*len(data_MV)] if prob[j][3] == 0 else ["--"]*len(data_MV),
                    ']': ciarkyMV, 'deltaMV_H': dHMV, '>': ciarkyMV, 
                    '1-norm E_MV,H': sol_4['x'][5*len(data_CV) + 3*len(data_MV):5*len(data_CV) + 4*len(data_MV)] if prob[j][3] == 0 else ["--"]*len(data_MV)}
            epsy4 = {'Variable': ['MV ' + str(data_MV[i]['MV Index']) for i in data_MV], 'Name': [data_MV[i]['MV Name'] for i in data_MV],
                    '||': ciarkyMV2, '-1-norm E_MV,L': -sol_4['x'][5*len(data_CV) + 4*len(data_MV):5*len(data_CV) + 5*len(data_MV)] if prob[j][3] == 0 else ["--"]*len(data_MV),
                    '<': ciarkyMV, '0': [0]*len_MV, '[': ciarkyMV,
                    '1-norm e_MV,L': sol_4['x'][3*len(data_CV) + 2*len(data_MV):3*len(data_CV) + 3*len(data_MV)] if prob[j][3] == 0 else ["--"]*len(data_MV),
                    ']': ciarkyMV, 'deltaMV_L': dLMV, '>': ciarkyMV,
                    '1-norm E_MV,L': sol_4['x'][5*len(data_CV) + 4*len(data_MV):5*len(data_CV) + 5*len(data_MV)] if prob[j][3] == 0 else ["--"]*len(data_MV)}

            # Display 1-norm epsilon values - .txt file
            with open("epsilons_1_norm_" + fileURT.split('.')[0] + model_type[j] + ".txt", "w") as f:
                f.write('\n=====================\n1-norm epsilon VALUES\n====================\n\n')
                f.write('\n(0 <= e_CV_H & -E_CV_H <= e_CV_H <= E_CV_H)\n\n')
                f.write(tabulate(pd.DataFrame(data = epsy1), headers = 'keys', showindex = False, floatfmt = ".4E"))
                f.write('\n\n\n(0 <= e_CV_H & -E_CV_H <= e_CV_H <= E_CV_H)\n\n')
                f.write(tabulate(pd.DataFrame(data = epsy2), headers = 'keys', showindex = False, floatfmt = ".4E"))
                f.write('\n\n\n(0 <= e_MV_H <= deltaMV_H & -E_MV_H <= e_MV_H <= E_MV_H)\n\n')
                f.write(tabulate(pd.DataFrame(data = epsy3), headers = 'keys', showindex = False, floatfmt = ".4E"))
                f.write('\n\n\n(0 <= e_MV_L <= deltaMV_L & -E_MV_L <= e_MV_L <= E_MV_L)\n\n')
                f.write(tabulate(pd.DataFrame(data = epsy4), headers = 'keys', showindex = False, floatfmt = ".4E"))
            if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Solve is True:
                os.system('start "" /max ' + 'epsilons_1_norm_' + fileURT.split('.')[0] + model_type[j] + '.txt')
            if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                print(), print("\nEpsilons - 1-norm\n==================\n")
                print(), print("(0 <= e_CV_H & -E_CV_H <= e_CV_H <= E_CV_H)\n")
                print(tabulate(epsy1, headers = 'keys', floatfmt = ".4E"))
                print(), print("(0 <= e_CV_H & -E_CV_H <= e_CV_H <= E_CV_H)\n")
                print(tabulate(epsy2, headers = 'keys', floatfmt = ".4E"))
                print(), print("(0 <= e_MV_H <= ΔMV_H & -E_MV_H <= e_MV_H <= E_MV_H)\n")
                print(tabulate(epsy3, headers = 'keys', floatfmt = ".4E"))
                print(), print("(0 <= e_MV_L <= ΔMV_L & -E_MV_L <= e_MV_L <= E_MV_L)\n")
                print(tabulate(epsy4, headers = 'keys', floatfmt = ".4E"))
                print()

            udaje = {'Info': ['Gap', 'Relative gap', 'Primal objective', 'Dual objective', 'Primal infeasibility', 'Dual infeasibility', 'Primal slack', 'Dual slack'],
                    'Hard bounds': [str(sol3['gap']) if prob[j][0] == 0 else "--", str(sol3['relative gap']) if prob[j][0] == 0 else "--", str(sol3['primal objective']) if prob[j][0] == 0 else "--", \
                                    str(sol3['dual objective']) if prob[j][0] == 0 else "--", str(sol3['primal infeasibility']) if prob[j][0] == 0 else "--", str(sol3['dual infeasibility']) if prob[j][0] == 0 else "--", \
                                    str(sol3['primal slack']) if prob[j][0] == 0 else "--", str(sol3['dual slack']) if prob[j][0] == 0 else "--"],
                    'Tight hard bounds': [str(sol['gap']) if prob[j][1] == 0 else "--", str(sol['relative gap']) if prob[j][1] == 0 else "--", str(sol['primal objective']) if prob[j][1] == 0 else "--", \
                                            str(sol['dual objective']) if prob[j][1] == 0 else "--", str(sol['primal infeasibility']) if prob[j][1] == 0 else "--", str(sol['dual infeasibility']) if prob[j][1] == 0 else "--", \
                                            str(sol['primal slack']) if prob[j][1] == 0 else "--", str(sol['dual slack']) if prob[j][1] == 0 else "--"],
                    '2-norm^2 soft bounds': [str(sol_2['gap']) if prob[j][2] == 0 else "--", str(sol_2['relative gap']) if prob[j][2] == 0 else "--", str(sol_2['primal objective']) if prob[j][2] == 0 else "--", \
                                                str(sol_2['dual objective']) if prob[j][2] == 0 else "--", str(sol_2['primal infeasibility']) if prob[j][2] == 0 else "--", str(sol_2['dual infeasibility']) if prob[j][2] == 0 else "--", \
                                                str(sol_2['primal slack']) if prob[j][2] == 0 else "--", str(sol_2['dual slack']) if prob[j][2] == 0 else "--"],
                    '1-norm soft bounds': [str(sol_4['gap']) if prob[j][3] == 0 else "--", str(sol_4['relative gap']) if prob[j][3] == 0 else "--", str(sol_4['primal objective']) if prob[j][3] == 0 else "--", \
                                            str(sol_4['dual objective']) if prob[j][3] == 0 else "--", str(sol_4['primal infeasibility']) if prob[j][3] == 0 else "--", str(sol_4['dual infeasibility']) if prob[j][3] == 0 else "--", \
                                            str(sol_4['primal slack']) if prob[j][3] == 0 else "--", str(sol_4['dual slack']) if prob[j][3] == 0 else "--"]}
            udaje2 = {'': ['Hard bounds', 'Tight hard bounds', '2-norm^2 soft bounds', '1-norm soft bounds'],
                    'Status': [sol_3['status'] if prob[j][0] == 0 else "--", sol['status'] if prob[j][1] == 0 else "--", sol_2['status'] if prob[j][2] == 0 else "--", sol_4['status'] if prob[j][3] == 0 else "--"]}
            with open("solution_data_" + fileURT.split('.')[0] + model_type[j] + ".txt", "w") as f:
                f.write('\n=====================\nSolution data\n====================\n\n')
                f.write(tabulate(pd.DataFrame(data = udaje2), headers = 'keys', showindex = False))
                f.write('\n\n')
                f.write(tabulate(pd.DataFrame(data = udaje), headers = 'keys', showindex = False))
            if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Solve is True:
                os.system('start "" /max ' + 'solution_data_' + fileURT.split('.')[0] + model_type[j] + '.txt')
            if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                print()
                print("\nSolution data\n=============\n")
                print(tabulate(udaje2, headers = 'keys'))
                print()
                print(tabulate(udaje, headers = 'keys'))
                print()

            if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                # Calculation time
                cas = {'HARD': [elapsed_BEZ_SOFTU],
                    'TIGHT HARD': [elapsed_PUK_POVODNE],
                    '2-norm^2': [elapsed_2_NORMA],
                    '1-norm': [elapsed_1_NORMA]};
                print("\nElapsed Time\n============\n")
                print(tabulate(cas, headers = 'keys', floatfmt = ".8f"))
                print()

        # Check solution status
        statuses = {'HARD BOUNDS real': sol3, 'TIGHT HARD BOUNDS real': sol, '2-norm^2 BOUNDS real': sol2, '1-norm BOUNDS real': sol4, 
                    'HARD BOUNDS ideal': sol_ideal3, 'TIGHT HARD BOUNDS ideal': sol_ideal, '2-norm^2 BOUNDS ideal': sol_ideal2, '1-norm BOUNDS ideal': sol_ideal4}
        descs = ['_tight_hard', '_2_norm', '_hard', '_1_norm', '_tight_hard', '_2_norm', '_hard', '_1_norm']
        for desc, stat in enumerate(statuses.keys()):
            if prob[int(desc/4)][desc%4] == 0:
                if statuses[stat]['status'] != 'optimal' and (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                    if 'BFCCU2C1_data' in sys.argv[-2]:
                        print("\n!!!! Solver did not find optimal solution for " + stat + " model !!!!\n\n")
                    else:
                        print('\n' + stylize("Solver did not find optimal solution for " + stat + " model", fg('dark_orange') + attr('bold')) + '\n\n')
                    table_name = 'Recommendations_' + fileURT.split('.')[0] + descs[desc] + '.txt'
                    with open(table_name, "w") as file:
                        message = '\nThere is no optimal solution. Gap: ' + str(statuses[stat]['gap']) + '\n\n'
                        file.write(message)
                        file.close()
                    if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Solve is True:
                        os.system('start "" /max ' + table_name)
            elif args_data.Solve is True:
                print("\n!!!! Solver converged to INF for " + stat + " model !!!!\n")
                table_name = 'Recommendations_' + fileURT.split('.')[0] + descs[desc] + '.txt'
                with open(table_name, "w") as file:
                    message = '\nThere is no optimal solution. Gap: INF\n\n'
                    file.write(message)
                    file.close()
                if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Solve is True:
                    os.system('start "" /max ' + table_name)

        # Solutions dicitonaries
        solution = {'real': sol, 'ideal': sol_ideal}
        solution2 = {'real': sol2, 'ideal': sol_ideal2}
        solution3 = {'real': sol3, 'ideal': sol_ideal3}
        solution4 = {'real': sol4, 'ideal': sol_ideal4}

        # Insert steady state data
        data_CV, data_MV = insert_solution_into_data(solution, solution2, solution3, solution4, data_CV, data_MV, prob, 2)

        # Calculation solved
        if args_data.Solve is True:
            if 'BFCCU2C1_data' in sys.argv[-2]:
                print('\n!!!! Calculation of optimal point solved !!!!\n\n')
            else:
                print('\n' + stylize('Calculation of optimal point solved', fg('green') + attr('bold')) + '\n\n')

        # Save data to external file
        if os.path.isfile("data_" + fileURT + "_" + fileExcel + "_S.txt") is True:
            t1, t2 = load_from_file("data_" + fileURT + "_" + fileExcel + "_S.txt", 0)
            (_, _, _, _, _, _, _, _, mtime_URT, _) = os.stat(fileURT)
            (_, _, _, _, _, _, _, _, mtime_XLSX, _) = os.stat(fileExcel)
            if t1 < int(mtime_URT) or t2 < int(mtime_XLSX):
                save_to_file(solution, solution2, solution3, solution4, M, M2, M3, M4, data_CV, data_MV, data, fileExcel, fileURT)
                if args_data.Solve is True:
                    if 'BFCCU2C1_data' in sys.argv[-2]:
                        print('\n!!!! Data was saved to external file !!!!\n\n')
                    else:
                        print('\n' + stylize('Data was saved to external file', fg('blue') + attr('bold')) + '\n\n')
        else:
            save_to_file(solution, solution2, solution3, solution4, M, M2, M3, M4, data_CV, data_MV, data, fileExcel, fileURT)
            if args_data.Solve is True:
                if 'BFCCU2C1_data' in sys.argv[-2]:
                    print('\n!!!! Data was saved to external file !!!!\n\n')
                else:
                    print('\n' + stylize('Data was saved to external file', fg('blue') + attr('bold')) + '\n\n')

        if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Solve is True:
            os.system('start "" /max ' + 'limits_' + fileURT.split('.')[0] + '_real.txt')
            os.system('start "" /max ' + 'limits_' + fileURT.split('.')[0] + '_ideal.txt')
            os.system('start "" /max ' + 'objective_' + fileURT.split('.')[0] + '_real.txt')
            os.system('start "" /max ' + 'objective_' + fileURT.split('.')[0] + '_ideal.txt')
            os.system('start "" /max ' + 'epsilons_2_norm_' + fileURT.split('.')[0] + '_real.txt')
            os.system('start "" /max ' + 'epsilons_2_norm_' + fileURT.split('.')[0] + '_ideal.txt')
            os.system('start "" /max ' + 'epsilons_1_norm_' + fileURT.split('.')[0] + '_real.txt')
            os.system('start "" /max ' + 'epsilons_1_norm_' + fileURT.split('.')[0] + '_ideal.txt')
            os.system('start "" /max ' + 'solution_data_' + fileURT.split('.')[0] + '_real.txt')
            os.system('start "" /max ' + 'solution_data_' + fileURT.split('.')[0] + '_ideal.txt')

    except Exception as e:
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)
        if 'BFCCU2C1_data' in sys.argv[-2]:
            print('\n!!!! Load correct files and solve problem !!!!\n\n')
        else:
            print('\n' + stylize('Load correct files and solve problem', fg('red') + attr('bold')) + '\n\n')

def execute2(args_data):

    # Create global variables
    global solution, solution2, solution3, solution4, M, M2, M3, M4, data_CV, data_MV, data
    M, M2, M3, M4 = {}, {}, {}, {}
    solution, solution2, solution3, solution4 = {}, {}, {}, {}

    try:
        # Load files
        fileExcel = re.search('Zošit(.+?).xlsx', args_data.xlsx).group(0)
        fileURT = re.search('BFCCU2C1_data(.+?).urt', args_data.urt).group(0)

        # Import data from Excel
        config_values = pd.read_excel(fileExcel, 'Sheet1')
        SS_values = pd.read_excel(fileExcel, 'Sheet2')

        # Import data from URT
        tree = ET.parse(fileURT)
        rootET = tree.getroot()

        # Load variables
        variables, place, gains, properties, _ = var_names()

        # Collect CVs, MVs & DVs data
        data = collect_data2(config_values, SS_values, args_data, rootET, place)

        # Collect gain matrix data
        data_gain_matrix = collect_data_for_gain_matrix(rootET, 'CV', 'Array of CV info structures', gains)

        # Load Model Usage data
        model_usage_vector_data = model_usage_vector(rootET)

        # Model gain position template
        model_usage_array = create_model_usage_array(data_gain_matrix, model_usage_vector_data, list(data.values())[0][1], args_data)

        # Load Model Usage data
        model_usage_vector_data = model_usage_vector(rootET)

        # Create gain matrix
        gain_matrix_model = create_gain_matrix(data_gain_matrix, gains)

        # Filter data by status
        data_CV, index_CV, data_MV, index_MV, data_DV, index_DV, gain_MV, gain_DV = {}, {}, {}, {}, {}, {}, {}, {}
        for num_X, dataX in data.items():
            data_CV[num_X], index_CV[num_X], data_MV[num_X], index_MV[num_X], data_DV[num_X], index_DV[num_X], data[num_X] = filter_data_by_status(dataX)

            # Adjust the gain matrix to the final form
            gain_MV[num_X], gain_DV[num_X] = filter_gain_matrix(index_CV[num_X], index_MV[num_X], index_DV[num_X], gain_matrix_model, model_usage_array)

        # Print saving data to file
        if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Solve is True:
            print(a[10])

        for dtX in range(1, len(data.items()) + 1):

            # Set file open method
            h_var = "w" if dtX == 1 else "a"

            if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                print("\n" + a[5] + "\n")
                print(config_values.columns[dtX].strftime("%Y-%m-%d %H:%M:%S"))
                print("\n" + a[5] + "\n")

            # Column separators
            len_CV, len_MV = len(data_CV[dtX]), len(data_MV[dtX])
            ciarkyCV, ciarkyMV, ciarkyCV2, ciarkyMV2, ciarkyCV3, ciarkyMV3 = ["|"]*len_CV, ["|"]*len_MV, ["||"]*len_CV, ["||"]*len_MV, ["|||"]*len_CV, ["|||"]*len_MV

            # Real & ideal high & low limits
            HCV = [data_CV[dtX][i]['High Limit'] for i in data_CV[dtX]]
            LCV = [data_CV[dtX][i]['Low Limit'] for i in data_CV[dtX]]
            HMV = [data_MV[dtX][i]['High Limit'] for i in data_MV[dtX]]
            LMV = [data_MV[dtX][i]['Low Limit'] for i in data_MV[dtX]]
            iHCV = [data_CV[dtX][i]['Ideal High'] for i in data_CV[dtX]]
            iLCV = [data_CV[dtX][i]['Ideal Low'] for i in data_CV[dtX]]
            iHMV = [data_MV[dtX][i]['Ideal High'] for i in data_MV[dtX]]
            iLMV = [data_MV[dtX][i]['Ideal Low'] for i in data_MV[dtX]]
            HLni = [[HCV, LCV, HMV, LMV], [iHCV, iLCV, iHMV, iLMV]]

            # Delta soft high & low limits
            dHCV = [data_CV[dtX][i]['Delta Soft High Limit'] for i in data_CV[dtX]]
            dLCV = [data_CV[dtX][i]['Delta Soft Low Limit'] for i in data_CV[dtX]]
            dHMV = [data_MV[dtX][i]['Delta Soft High Limit'] for i in data_MV[dtX]]
            dLMV = [data_MV[dtX][i]['Delta Soft Low Limit'] for i in data_MV[dtX]]
            dHLni = [[(np.array(HCV) - np.array(dHCV)).tolist(), (np.array(LCV) + np.array(dLCV)).tolist(), (np.array(HMV) - np.array(dHMV)).tolist(), \
                    (np.array(LMV) + np.array(dLMV)).tolist()], [(np.array(iHCV) - np.array(dHCV)).tolist(), (np.array(iLCV) + np.array(dLCV)).tolist(), \
                    (np.array(iHMV) - np.array(dHMV)).tolist(), (np.array(iLMV) + np.array(dLMV)).tolist()]]

            # Matrices (P, q, r, G, h - ideal & real, A, b)
            P, P2, P4 = create_matrix_P(data_CV[dtX], data_MV[dtX])
            q, q2, q4 = create_matrix_q(data_CV[dtX], data_MV[dtX])
            G, G2, G4 = create_matrix_G(data_CV[dtX], data_MV[dtX])
            h, h2, h3, h4, data_CV[dtX], data_MV[dtX] = create_matrix_h(data_CV[dtX], data_MV[dtX])
            h_ideal, h_ideal2, h_ideal3, h_ideal4, data_CV[dtX], data_MV[dtX] = create_matrix_h(data_CV[dtX], data_MV[dtX], ideal = 1)
            A, A2, A4 = create_matrix_A(data_CV[dtX], gain_MV[dtX], data_MV[dtX])
            b = create_matrix_b(data_CV[dtX], data_MV[dtX], gain_MV[dtX])
            r = create_matrix_r(data_CV[dtX], data_MV[dtX])

            #   - PUK
            # 2 - 2-norm^2
            # 3 - HARD
            # 4 - 1-norm

            # Matrices dictionaries
            M[dtX] = {'P': P, 'q': q, 'G': G, 'h': h, 'A': A, 'b': b, 'r': r, 'h_ideal': h_ideal}
            M2[dtX] = {'P': P2, 'q': q2, 'G': G2, 'h': h2, 'A': A2, 'b': b, 'r': r, 'h_ideal': h_ideal2}
            M3[dtX] = {'P': P, 'q': q, 'G': G, 'h': h3, 'A': A, 'b': b, 'r': r, 'h_ideal': h_ideal3}
            M4[dtX] = {'P': P4, 'q': q4, 'G': G4, 'h': h4, 'A': A4, 'b': b, 'r': r, 'h_ideal': h_ideal4}

            # Solve optimization problems (optimal steady states) - real models
            prob_R = [0,0,0,0]
            t_BEZ_SOFTU = time.time();
            try:
                sol3 = solvers.qp(2*P, q, G, h3, A, b, solver = 'mosek')
                if sol3['status'] != "optimal":
                    sol3 = 0
                    prob_R[0] = 1
            except:
                sol3 = 0
                prob_R[0] = 1
            elapsed_BEZ_SOFTU = time.time() - t_BEZ_SOFTU;
            t_PUK_POVODNE = time.time();
            try:
                sol = solvers.qp(2*P, q, G, h, A, b, solver = 'mosek')
                if sol['status'] != "optimal":
                    sol = 0
                    prob_R[1] = 1
            except:
                sol = 0
                prob_R[1] = 1
            elapsed_PUK_POVODNE = time.time() - t_PUK_POVODNE;
            t_2_NORMA = time.time();
            try:
                sol2 = solvers.qp(2*P2, q2, G2, h2, A2, b, solver = 'mosek')
                if sol2['status'] != "optimal":
                    sol2 = 0
                    prob_R[2] = 1
            except:
                sol2 = 0
                prob_R[2] = 1
            elapsed_2_NORMA = time.time() - t_2_NORMA;
            t_1_NORMA = time.time();
            try:
                sol4 = solvers.qp(2*P4, q4, G4, h4, A4, b, solver = 'mosek')
                if sol4['status'] != "optimal":
                    sol4 = 0
                    prob_R[3] = 1
            except:
                sol4 = 0
                prob_R[3] = 1
            elapsed_1_NORMA = time.time() - t_1_NORMA;

            # Solve optimization problems (optimal steady states) - ideal models
            prob_I = [0,0,0,0]
            try:
                sol_ideal3 = solvers.qp(2*P, q, G, h_ideal3, A, b, solver = 'mosek')
                if sol_ideal3['status'] != "optimal":
                    sol_ideal3 = 0
                    prob_I[0] = 1
            except:
                sol_ideal3 = 0
                prob_I[0] = 1
            try:
                sol_ideal = solvers.qp(2*P, q, G, h_ideal, A, b, solver = 'mosek')
                if sol_ideal['status'] != "optimal":
                    sol_ideal = 0
                    prob_I[1] = 1
            except:
                sol_ideal = 0
                prob_I[1] = 1
            try:
                sol_ideal2 = solvers.qp(2*P2, q2, G2, h_ideal2, A2, b, solver = 'mosek')
                if sol_ideal2['status'] != "optimal":
                    sol_ideal2 = 0
                    prob_I[2] = 1
            except:
                sol_ideal2 = 0
                prob_I[2] = 1
            try:
                sol_ideal4 = solvers.qp(2*P4, q4, G4, h_ideal4, A4, b, solver = 'mosek')
                if sol_ideal4['status'] != "optimal":
                    sol_ideal4 = 0
                    prob_I[3] = 1
            except:
                sol_ideal4 = 0
                prob_I[3] = 1

            print_type, model_type, prob = [a[27], a[28]], ['_real', '_ideal'], [prob_R, prob_I]
            model_data = [[sol, sol2, sol3, sol4], [sol_ideal, sol_ideal2, sol_ideal3, sol_ideal4]]

            for j, i in enumerate(model_data):
                sol_, sol_2, sol_3, sol_4 = i[0], i[1], i[2], i[3]
                h_var = "w" if dtX == 1 else "a"

                # Model type
                if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                    print(print_type[j])

                # Limit violations
                if prob_R[0] == 0 or prob_I[0] == 0:
                    CV_L_S1 = (np.array([round(float(i - j), 2) > 0.0 for i, j in zip(np.array(dHLni[j][1]), np.array(sol_3['x'][:len(data_CV[dtX])]))]) > 0).tolist()
                    CV_H_S1 = (np.array([round(float(i - j), 2) < 0.0 for i, j in zip(np.array(dHLni[j][0]), np.array(sol_3['x'][:len(data_CV[dtX])]))]) > 0).tolist()
                    MV_L_S1 = (np.array([round(float(i - j), 2) > 0.0 for i, j in zip(np.array(dHLni[j][3]), np.array(sol_3['x'][len(data_CV[dtX]):len(data_CV[dtX]) + len(data_MV[dtX])]))]) > 0).tolist()
                    MV_H_S1 = (np.array([round(float(i - j), 2) < 0.0 for i, j in zip(np.array(dHLni[j][2]), np.array(sol_3['x'][len(data_CV[dtX]):len(data_CV[dtX]) + len(data_MV[dtX])]))]) > 0).tolist()
                else:
                    CV_L_S1, CV_H_S1, MV_L_S1, MV_H_S1 = len(data_CV[dtX])*["--"], len(data_CV[dtX])*["--"], len(data_MV[dtX])*["--"], len(data_MV[dtX])*["--"]
                if prob_R[1] == 0 or prob_I[1] == 0:
                    CV_L_S2 = (np.array([round(float(i - j), 2) > 0.0 for i, j in zip(np.array(dHLni[j][1]), np.array(sol_['x'][:len(data_CV[dtX])]))]) > 0).tolist()
                    CV_H_S2 = (np.array([round(float(i - j), 2) < 0.0 for i, j in zip(np.array(dHLni[j][0]), np.array(sol_['x'][:len(data_CV[dtX])]))]) > 0).tolist()
                    MV_L_S2 = (np.array([round(float(i - j), 2) > 0.0 for i, j in zip(np.array(dHLni[j][3]), np.array(sol_['x'][len(data_CV[dtX]):len(data_CV[dtX]) + len(data_MV[dtX])]))]) > 0).tolist()
                    MV_H_S2 = (np.array([round(float(i - j), 2) < 0.0 for i, j in zip(np.array(dHLni[j][2]), np.array(sol_['x'][len(data_CV[dtX]):len(data_CV[dtX]) + len(data_MV[dtX])]))]) > 0).tolist()
                else:
                    CV_L_S2, CV_H_S2, MV_L_S2, MV_H_S2 = len(data_CV[dtX])*["--"], len(data_CV[dtX])*["--"], len(data_MV[dtX])*["--"], len(data_MV[dtX])*["--"]
                if prob_R[2] == 0 or prob_I[2] == 0:
                    CV_L_S3 = (np.array([round(float(i - j), 2) > 0.0 for i, j in zip(np.array(dHLni[j][1]), np.array(sol_2['x'][:len(data_CV[dtX])]))]) > 0).tolist()
                    CV_H_S3 = (np.array([round(float(i - j), 2) < 0.0 for i, j in zip(np.array(dHLni[j][0]), np.array(sol_2['x'][:len(data_CV[dtX])]))]) > 0).tolist()
                    MV_L_S3 = (np.array([round(float(i - j), 2) > 0.0 for i, j in zip(np.array(dHLni[j][3]), np.array(sol_2['x'][len(data_CV[dtX]):len(data_CV[dtX]) + len(data_MV[dtX])]))]) > 0).tolist()
                    MV_H_S3 = (np.array([round(float(i - j), 2) < 0.0 for i, j in zip(np.array(dHLni[j][2]), np.array(sol_2['x'][len(data_CV[dtX]):len(data_CV[dtX]) + len(data_MV[dtX])]))]) > 0).tolist()
                else:
                    CV_L_S3, CV_H_S3, MV_L_S3, MV_H_S3 = len(data_CV[dtX])*["--"], len(data_CV[dtX])*["--"], len(data_MV[dtX])*["--"], len(data_MV[dtX])*["--"]
                if prob_R[3] == 0 or prob_I[3] == 0:
                    CV_L_S4 = (np.array([round(float(i - j), 2) > 0.0 for i, j in zip(np.array(dHLni[j][1]), np.array(sol_4['x'][:len(data_CV[dtX])]))]) > 0).tolist()
                    CV_H_S4 = (np.array([round(float(i - j), 2) < 0.0 for i, j in zip(np.array(dHLni[j][0]), np.array(sol_4['x'][:len(data_CV[dtX])]))]) > 0).tolist()
                    MV_L_S4 = (np.array([round(float(i - j), 2) > 0.0 for i, j in zip(np.array(dHLni[j][3]), np.array(sol_4['x'][len(data_CV[dtX]):len(data_CV[dtX]) + len(data_MV[dtX])]))]) > 0).tolist()
                    MV_H_S4 = (np.array([round(float(i - j), 2) < 0.0 for i, j in zip(np.array(dHLni[j][2]), np.array(sol_4['x'][len(data_CV[dtX]):len(data_CV[dtX]) + len(data_MV[dtX])]))]) > 0).tolist()
                else:
                    CV_L_S4, CV_H_S4, MV_L_S4, MV_H_S4 = len(data_CV[dtX])*["--"], len(data_CV[dtX])*["--"], len(data_MV[dtX])*["--"], len(data_MV[dtX])*["--"]

                CV_1, CV_2, CV_3, CV_4, MV_1, MV_2, MV_3, MV_4 = [], [], [], [], [], [], [], []
                for i in range(len(CV_L_S1)):
                    CV_1.append("LL/HL") if CV_L_S1[i] is True and CV_H_S1[i] is True else CV_1.append("LL") if CV_L_S1[i] is True else CV_1.append("HL") if CV_H_S1[i] is True else CV_1.append("--")
                    CV_2.append("LL/HL") if CV_L_S2[i] is True and CV_H_S2[i] is True else CV_2.append("LL") if CV_L_S2[i] is True else CV_2.append("HL") if CV_H_S2[i] is True else CV_2.append("--")
                    CV_3.append("LL/HL") if CV_L_S3[i] is True and CV_H_S3[i] is True else CV_3.append("LL") if CV_L_S3[i] is True else CV_3.append("HL") if CV_H_S3[i] is True else CV_3.append("--")
                    CV_4.append("LL/HL") if CV_L_S4[i] is True and CV_H_S4[i] is True else CV_4.append("LL") if CV_L_S4[i] is True else CV_4.append("HL") if CV_H_S4[i] is True else CV_4.append("--")
                for i in range(len(MV_L_S1)):
                    MV_1.append("LL/HL") if MV_L_S1[i] is True and MV_H_S1[i] is True else MV_1.append("LL") if MV_L_S1[i] is True else MV_1.append("HL") if MV_H_S1[i] is True else MV_1.append("--")
                    MV_2.append("LL/HL") if MV_L_S2[i] is True and MV_H_S2[i] is True else MV_2.append("LL") if MV_L_S2[i] is True else MV_2.append("HL") if MV_H_S2[i] is True else MV_2.append("--")
                    MV_3.append("LL/HL") if MV_L_S3[i] is True and MV_H_S3[i] is True else MV_3.append("LL") if MV_L_S3[i] is True else MV_3.append("HL") if MV_H_S3[i] is True else MV_3.append("--")
                    MV_4.append("LL/HL") if MV_L_S4[i] is True and MV_H_S4[i] is True else MV_4.append("LL") if MV_L_S4[i] is True else MV_4.append("HL") if MV_H_S4[i] is True else MV_4.append("--")                

                # Display limit data - command line
                cv_val = {'Variable': ['CV ' + str(data_CV[dtX][i]['CV Index']) for i in data_CV[dtX]], 'Name': [data_CV[dtX][i]['CV Name'] for i in data_CV[dtX]],
                        '||': ciarkyCV2, 'Low hard limit': HLni[j][1], 'Low soft limit': dHLni[j][1], '[': ciarkyCV,
                        'Hard': sol_3['x'][:len(data_CV[dtX])] if prob[j][0] == 0 else ["--"]*len(data_CV[dtX]),
                        'Tight hard': sol_['x'][:len(data_CV[dtX])] if prob[j][1] == 0 else ["--"]*len(data_CV[dtX]),
                        '2-norm^2': sol_2['x'][:len(data_CV[dtX])] if prob[j][2] == 0 else ["--"]*len(data_CV[dtX]),
                        '1-norm': sol_4['x'][:len(data_CV[dtX])] if prob[j][3] == 0 else ["--"]*len(data_CV[dtX]),
                        ']': ciarkyCV, 'High soft limit': dHLni[j][0], 'High hard limit': HLni[j][0], '|||': ciarkyCV3,
                        'VIO Hard': CV_1, 'VIO Tight hard': CV_2, 'VIO 2-norm^2': CV_3, 'VIO 1-norm': CV_4}
                mv_val = {'Variable': ['MV ' + str(data_MV[dtX][i]['MV Index']) for i in data_MV[dtX]], 'Name': [data_MV[dtX][i]['MV Name'] for i in data_MV[dtX]],
                        '||': ciarkyMV2, 'Low hard limit': HLni[j][3], 'Low soft limit': dHLni[j][3], '[': ciarkyMV,
                        'Hard': sol_3['x'][len(data_CV[dtX]):len(data_CV[dtX]) + len(data_MV[dtX])] if prob[j][0] == 0 else ["--"]*len(data_MV[dtX]),
                        'Tight hard': sol_['x'][len(data_CV[dtX]):len(data_CV[dtX]) + len(data_MV[dtX])] if prob[j][1] == 0 else ["--"]*len(data_MV[dtX]),
                        '2-norm^2': sol_2['x'][len(data_CV[dtX]):len(data_CV[dtX]) + len(data_MV[dtX])] if prob[j][2] == 0 else ["--"]*len(data_MV[dtX]),
                        '1-norm': sol_4['x'][len(data_CV[dtX]):len(data_CV[dtX]) + len(data_MV[dtX])] if prob[j][3] == 0 else ["--"]*len(data_MV[dtX]),
                        ']': ciarkyMV, 'High soft limit': dHLni[j][2], 'High hard limit': HLni[j][2], '|||': ciarkyMV3,
                        'VIO Hard': MV_1, 'VIO Tight hard': MV_2, 'VIO 2-norm^2': MV_3, 'VIO 1-norm': MV_4}

                # Display limit data - .txt file
                with open("limits_" + fileURT.split('.')[0] + model_type[j] + "_M.txt", h_var) as f:
                    if dtX == 1: f.write(a[21]), f.write(a[22])
                    if dtX != 1: f.write('\n')
                    f.write('\n\n' + a[4] + str(config_values.columns[dtX].strftime("%Y-%m-%d %H:%M:%S")) + a[4])
                    f.write('\n==============\nCV constraints\n==============\n\n')
                    f.write(tabulate(pd.DataFrame(data = cv_val), headers = 'keys', showindex = False, floatfmt = ".2f"))
                    f.write('\n\n\n==============\nMV constraints\n==============\n\n')
                    f.write(tabulate(pd.DataFrame(data = mv_val), headers = 'keys', showindex = False, floatfmt = ".2f"))
                if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                    print()
                    print("\nCV optimal\n==========\n")
                    print(tabulate(cv_val, headers = 'keys', floatfmt = ".2f"))
                    print()
                    print("\nMV optimal\n==========\n")
                    print(tabulate(mv_val, headers = 'keys', floatfmt = ".2f"))
                    print()

                # Display objective function data
                obj = {'HARD': sol_3['x'].trans()*P*sol_3['x'] + q.trans()*sol_3['x'] + r if prob[j][0] == 0 else "----",
                    'TIGHT HARD': sol_['x'].trans()*P*sol_['x'] + q.trans()*sol_['x'] + r if prob[j][1] == 0 else "----",
                    '2-norm^2 WITH e': sol_2['x'].trans()*P2*sol_2['x'] + q2.trans()*sol_2['x'] + r if prob[j][2] == 0 else "----",
                    '2-norm^2 WITHOUT e': sol_2['x'][:len(data_CV[dtX]) + len(data_MV[dtX])].trans()*P2[:len(data_CV[dtX]) + len(data_MV[dtX]),:len(data_CV[dtX]) + len(data_MV[dtX])]*sol_2['x'][:len(data_CV[dtX]) + len(data_MV[dtX])] + q2[:len(data_CV[dtX]) + len(data_MV[dtX])].trans()*sol_2['x'][:len(data_CV[dtX]) + len(data_MV[dtX])] + r if prob[j][2] == 0 else "----",
                    '1-norm WITH e': sol_4['x'].trans()*P4*sol_4['x'] + q4.trans()*sol_4['x'] + r if prob[j][3] == 0 else "----",
                    '1-norm WITHOUT e': sol_4['x'][:len(data_CV[dtX]) + len(data_MV[dtX])].trans()*P4[:len(data_CV[dtX]) + len(data_MV[dtX]),:len(data_CV[dtX]) + len(data_MV[dtX])]*sol_4['x'][:len(data_CV[dtX]) + len(data_MV[dtX])] + q4[:len(data_CV[dtX]) + len(data_MV[dtX])].trans()*sol_4['x'][:len(data_CV[dtX]) + len(data_MV[dtX])] + r if prob[j][3] == 0 else "----"}
                if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                    print()
                    print("\nObjective function value\n========================\n")
                    print(tabulate(obj, headers = 'keys', floatfmt = ".4E"))
                    print()

                # Display relative objective function data
                rozdiel = {'HARD': (obj['HARD'] - obj['HARD'])/abs(obj['HARD'])*100 if str(obj['HARD']) != "----" else matrix(np.array([-1.0])),
                        'TIGHT HARD': (obj['TIGHT HARD'] - obj['HARD'])/abs(obj['HARD'])*100 if (str(obj['HARD']) != "----" and str(obj['TIGHT HARD']) != "----") else matrix(np.array([-1.0])),
                        '2-norm^2 WITH e': (obj['2-norm^2 WITH e'] - obj['HARD'])/abs(obj['HARD'])*100 if (str(obj['HARD']) != "----" and str(obj['2-norm^2 WITH e']) != "----") else matrix(np.array([-1.0])),
                        '2-norm^2 WITHOUT e': (obj['2-norm^2 WITHOUT e'] - obj['HARD'])/abs(obj['HARD'])*100 if (str(obj['HARD']) != "----" and str(obj['2-norm^2 WITHOUT e']) != "----") else matrix(np.array([-1.0])),
                        '1-norm WITH e': (obj['1-norm WITH e'] - obj['HARD'])/abs(obj['HARD'])*100 if (str(obj['HARD']) != "----" and str(obj['1-norm WITH e']) != "----") else matrix(np.array([-1.0])),
                        '1-norm WITHOUT e': (obj['1-norm WITHOUT e'] - obj['HARD'])/abs(obj['HARD'])*100 if (str(obj['HARD']) != "----" and str(obj['1-norm WITHOUT e']) != "----") else matrix(np.array([-1.0]))}
                if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                    print()
                    print("\nRelative objective function value (against HARD in %)\n=====================================================\n")
                    print(tabulate(rozdiel, headers = 'keys', floatfmt = ".4E"))
                    print()

                # Display objective function data - .txt file
                with open("objective_" + fileURT.split('.')[0] + model_type[j] + "_M.txt", h_var) as f:
                    if dtX == 1: f.write(a[23]), f.write(a[24])
                    if dtX != 1: f.write('\n')
                    f.write('\n\n' + a[4] + str(config_values.columns[dtX].strftime("%Y-%m-%d %H:%M:%S")) + a[4])
                    f.write('\n===================\nObjective functions\n===================\n\n')
                    f.write(tabulate(pd.DataFrame(data = obj), headers = 'keys', showindex = False, floatfmt = ".4E"))
                    f.write('\n\n\n===================================\nRelative objective functions (in %)\n===================================\n\n')
                    f.write(tabulate(pd.DataFrame(data = rozdiel), headers = 'keys', showindex = False, floatfmt = ".4E"))

                # Display 2-norm^2 epsilon values
                epsy1 = {'Variable': ['CV ' + str(data_CV[dtX][i]['CV Index']) for i in data_CV[dtX]], 'Name': [data_CV[dtX][i]['CV Name'] for i in data_CV[dtX]],
                        '||': ciarkyCV2, '0': [0]*len_CV, '[': ciarkyCV,
                        'e_CV,H': sol_2['x'][len(data_CV[dtX]) + len(data_MV[dtX]):2*len(data_CV[dtX]) + len(data_MV[dtX])] if prob[j][2] == 0 else ["--"]*len(data_CV[dtX]),
                        ']': ciarkyCV, 'INF': ['inf']*len_CV}
                epsy2 = {'Variable': ['CV ' + str(data_CV[dtX][i]['CV Index']) for i in data_CV[dtX]], 'Name': [data_CV[dtX][i]['CV Name'] for i in data_CV[dtX]],
                        '||': ciarkyCV2, '0': [0]*len_CV, '[': ciarkyCV,
                        'e_CV,L': sol_2['x'][2*len(data_CV[dtX]) + len(data_MV[dtX]):3*len(data_CV[dtX]) + len(data_MV[dtX])] if prob[j][2] == 0 else ["--"]*len(data_CV[dtX]),
                        ']': ciarkyCV, 'INF': ['inf']*len_CV}
                epsy3 = {'Variable': ['MV ' + str(data_MV[dtX][i]['MV Index']) for i in data_MV[dtX]], 'Name': [data_MV[dtX][i]['MV Name'] for i in data_MV[dtX]],
                        '||': ciarkyMV2, '0': [0]*len_MV, '[': ciarkyMV,
                        'e_MV,H': sol_2['x'][3*len(data_CV[dtX]) + len(data_MV[dtX]):3*len(data_CV[dtX]) + 2*len(data_MV[dtX])] if prob[j][2] == 0 else ["--"]*len(data_MV[dtX]),
                        ']': ciarkyMV, 'deltaMV_H': dHMV}
                epsy4 = {'Variable': ['MV ' + str(data_MV[dtX][i]['MV Index']) for i in data_MV[dtX]], 'Name': [data_MV[dtX][i]['MV Name'] for i in data_MV[dtX]],
                        '||': ciarkyMV2, '0': [0]*len_MV, '[': ciarkyMV,
                        'e_MV,L': sol_2['x'][3*len(data_CV[dtX]) + 2*len(data_MV[dtX]):3*len(data_CV[dtX]) + 3*len(data_MV[dtX])] if prob[j][2] == 0 else ["--"]*len(data_MV[dtX]),
                        ']': ciarkyMV, 'deltaMV_L': dLMV}

                # Display 2-norm^2 epsilon values - .txt file
                with open("epsilons_2_norm_" + fileURT.split('.')[0] + model_type[j] + "_M.txt", h_var) as f:
                    if dtX == 1: f.write(a[25])
                    if dtX != 1: f.write('\n')
                    f.write('\n\n' + a[4] + str(config_values.columns[dtX].strftime("%Y-%m-%d %H:%M:%S")) + a[4])
                    f.write('\n=================\nEpsilons 2-norm^2\n=================\n\n')
                    f.write('\n(0 <= e_CV_H)\n\n')
                    f.write(tabulate(pd.DataFrame(data = epsy1), headers = 'keys', showindex = False, floatfmt = ".4E"))
                    f.write('\n\n\n(0 <= e_CV_L)\n\n')
                    f.write(tabulate(pd.DataFrame(data = epsy2), headers = 'keys', showindex = False, floatfmt = ".4E"))
                    f.write('\n\n\n(0 <= e_MV_H <= deltaMV_H)\n\n')
                    f.write(tabulate(pd.DataFrame(data = epsy3), headers = 'keys', showindex = False, floatfmt = ".4E"))
                    f.write('\n\n\n(0 <= e_MV_L <= deltaMV_L)\n\n')
                    f.write(tabulate(pd.DataFrame(data = epsy4), headers = 'keys', showindex = False, floatfmt = ".4E"))
                if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                    print(), print("\nEpsilons - 2-norm^2\n===================\n")
                    print(), print("(0 <= e_CV_H)\n")
                    print(tabulate(epsy1, headers = 'keys', floatfmt = ".4E"))
                    print(), print("(0 <= e_CV_L)\n")
                    print(tabulate(epsy2, headers = 'keys', floatfmt = ".4E"))
                    print(), print("(0 <= e_MV_H <= ΔMV_H)\n")
                    print(tabulate(epsy3, headers = 'keys', floatfmt = ".4E"))
                    print(), print("(0 <= e_MV_L <= ΔMV_L)\n")
                    print(tabulate(epsy4, headers = 'keys', floatfmt = ".4E"))
                    print()

                # Display 1-norm epsilon values - .txt file
                epsy1 = {'Variable': ['CV ' + str(data_CV[dtX][i]['CV Index']) for i in data_CV[dtX]], 'Name': [data_CV[dtX][i]['CV Name'] for i in data_CV[dtX]],
                        '||': ciarkyCV2, '-E_CV,H': sol_4['x'][3*len(data_CV[dtX]) + 3*len(data_MV[dtX]):4*len(data_CV[dtX]) + 3*len(data_MV[dtX])] if prob[j][3] == 0 else ["--"]*len(data_CV[dtX]),
                        '<': ciarkyCV, '0': [0]*len_CV, '[': ciarkyCV,
                        'e_CV,H': sol_4['x'][len(data_CV[dtX]) + len(data_MV[dtX]):2*len(data_CV[dtX]) + len(data_MV[dtX])] if prob[j][3] == 0 else ["--"]*len(data_CV[dtX]),
                        ']': ciarkyCV, 'INF': ['inf']*len_CV, '>': ciarkyCV,
                        'E_CV,H': -sol_4['x'][3*len(data_CV[dtX]) + 3*len(data_MV[dtX]):4*len(data_CV[dtX]) + 3*len(data_MV[dtX])] if prob[j][3] == 0 else ["--"]*len(data_CV[dtX])}
                epsy2 = {'Variable': ['CV ' + str(data_CV[dtX][i]['CV Index']) for i in data_CV[dtX]], 'Name': [data_CV[dtX][i]['CV Name'] for i in data_CV[dtX]],
                        '||': ciarkyCV2, '-E_CV,L': sol_4['x'][4*len(data_CV[dtX]) + 3*len(data_MV[dtX]):5*len(data_CV[dtX]) + 3*len(data_MV[dtX])] if prob[j][3] == 0 else ["--"]*len(data_CV[dtX]),
                        '<': ciarkyCV, '0': [0]*len_CV, '[': ciarkyCV,
                        'e_CV,L': sol_4['x'][2*len(data_CV[dtX]) + len(data_MV[dtX]):3*len(data_CV[dtX]) + len(data_MV[dtX])] if prob[j][3] == 0 else ["--"]*len(data_CV[dtX]),
                        ']': ciarkyCV, 'INF': ['inf']*len_CV, '>': ciarkyCV,
                        'E_CV,L': -sol_4['x'][4*len(data_CV[dtX]) + 3*len(data_MV[dtX]):5*len(data_CV[dtX]) + 3*len(data_MV[dtX])] if prob[j][3] == 0 else ["--"]*len(data_CV[dtX])}
                epsy3 = {'Variable': ['MV ' + str(data_MV[dtX][i]['MV Index']) for i in data_MV[dtX]], 'Name': [data_MV[dtX][i]['MV Name'] for i in data_MV[dtX]],
                        '||': ciarkyMV2, '-E_MV,H': -sol_4['x'][5*len(data_CV[dtX]) + 3*len(data_MV[dtX]):5*len(data_CV[dtX]) + 4*len(data_MV[dtX])] if prob[j][3] == 0 else ["--"]*len(data_MV[dtX]),
                        '<': ciarkyMV, '0': [0]*len_MV, '[': ciarkyMV,
                        'e_MV,H': sol_4['x'][3*len(data_CV[dtX]) + len(data_MV[dtX]):3*len(data_CV[dtX]) + 2*len(data_MV[dtX])] if prob[j][3] == 0 else ["--"]*len(data_MV[dtX]),
                        ']': ciarkyMV, 'deltaMV_H': dHMV, '>': ciarkyMV, 
                        'E_MV,H': sol_4['x'][5*len(data_CV[dtX]) + 3*len(data_MV[dtX]):5*len(data_CV[dtX]) + 4*len(data_MV[dtX])] if prob[j][3] == 0 else ["--"]*len(data_MV[dtX])}
                epsy4 = {'Variable': ['MV ' + str(data_MV[dtX][i]['MV Index']) for i in data_MV[dtX]], 'Name': [data_MV[dtX][i]['MV Name'] for i in data_MV[dtX]],
                        '||': ciarkyMV2, '-E_MV,L': -sol_4['x'][5*len(data_CV[dtX]) + 4*len(data_MV[dtX]):5*len(data_CV[dtX]) + 5*len(data_MV[dtX])] if prob[j][3] == 0 else ["--"]*len(data_MV[dtX]),
                        '<': ciarkyMV, '0': [0]*len_MV, '[': ciarkyMV,
                        'e_MV,L': sol_4['x'][3*len(data_CV[dtX]) + 2*len(data_MV[dtX]):3*len(data_CV[dtX]) + 3*len(data_MV[dtX])] if prob[j][3] == 0 else ["--"]*len(data_MV[dtX]),
                        ']': ciarkyMV, 'deltaMV_L': dLMV, '>': ciarkyMV,
                        'E_MV,L': sol_4['x'][5*len(data_CV[dtX]) + 4*len(data_MV[dtX]):5*len(data_CV[dtX]) + 5*len(data_MV[dtX])] if prob[j][3] == 0 else ["--"]*len(data_MV[dtX])}

                # Display 1-norm epsilon values - .txt file
                with open("epsilons_1_norm_" + fileURT.split('.')[0] + model_type[j] + "_M.txt", h_var) as f:
                    if dtX == 1: f.write(a[26])
                    if dtX != 1: f.write('\n')
                    f.write('\n\n' + a[4] + str(config_values.columns[dtX].strftime("%Y-%m-%d %H:%M:%S")) + a[4])
                    f.write('\n===============\nEpsilons 1-norm\n===============\n\n')
                    f.write('\n(0 <= e_CV_H & -E_CV_H <= e_CV_H <= E_CV_H)\n\n')
                    f.write(tabulate(pd.DataFrame(data = epsy1), headers = 'keys', showindex = False, floatfmt = ".4E"))
                    f.write('\n\n\n(0 <= e_CV_H & -E_CV_H <= e_CV_H <= E_CV_H)\n\n')
                    f.write(tabulate(pd.DataFrame(data = epsy2), headers = 'keys', showindex = False, floatfmt = ".4E"))
                    f.write('\n\n\n(0 <= e_MV_H <= deltaMV_H & -E_MV_H <= e_MV_H <= E_MV_H)\n\n')
                    f.write(tabulate(pd.DataFrame(data = epsy3), headers = 'keys', showindex = False, floatfmt = ".4E"))
                    f.write('\n\n\n(0 <= e_MV_L <= deltaMV_L & -E_MV_L <= e_MV_L <= E_MV_L)\n\n')
                    f.write(tabulate(pd.DataFrame(data = epsy4), headers = 'keys', showindex = False, floatfmt = ".4E"))
                if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                    print(), print("\nEpsilons - 1-norm\n==================\n")
                    print(), print("(0 <= e_CV_H & -E_CV_H <= e_CV_H <= E_CV_H)\n")
                    print(tabulate(epsy1, headers = 'keys', floatfmt = ".4E"))
                    print(), print("(0 <= e_CV_H & -E_CV_H <= e_CV_H <= E_CV_H)\n")
                    print(tabulate(epsy2, headers = 'keys', floatfmt = ".4E"))
                    print(), print("(0 <= e_MV_H <= ΔMV_H & -E_MV_H <= e_MV_H <= E_MV_H)\n")
                    print(tabulate(epsy3, headers = 'keys', floatfmt = ".4E"))
                    print(), print("(0 <= e_MV_L <= ΔMV_L & -E_MV_L <= e_MV_L <= E_MV_L)\n")
                    print(tabulate(epsy4, headers = 'keys', floatfmt = ".4E"))
                    print()

                udaje = {'Info': ['Gap', 'Relative gap', 'Primal objective', 'Dual objective', 'Primal infeasibility', 'Dual infeasibility', 'Primal slack', 'Dual slack'],
                        'Hard': [str(sol3['gap']) if prob[j][0] == 0 else "--", str(sol3['relative gap']) if prob[j][0] == 0 else "--", str(sol3['primal objective']) if prob[j][0] == 0 else "--", \
                                        str(sol3['dual objective']) if prob[j][0] == 0 else "--", str(sol3['primal infeasibility']) if prob[j][0] == 0 else "--", str(sol3['dual infeasibility']) if prob[j][0] == 0 else "--", \
                                        str(sol3['primal slack']) if prob[j][0] == 0 else "--", str(sol3['dual slack']) if prob[j][0] == 0 else "--"],
                        'Tight hard': [str(sol['gap']) if prob[j][1] == 0 else "--", str(sol['relative gap']) if prob[j][1] == 0 else "--", str(sol['primal objective']) if prob[j][1] == 0 else "--", \
                                              str(sol['dual objective']) if prob[j][1] == 0 else "--", str(sol['primal infeasibility']) if prob[j][1] == 0 else "--", str(sol['dual infeasibility']) if prob[j][1] == 0 else "--", \
                                              str(sol['primal slack']) if prob[j][1] == 0 else "--", str(sol['dual slack']) if prob[j][1] == 0 else "--"],
                        '2-norm^2': [str(sol_2['gap']) if prob[j][2] == 0 else "--", str(sol_2['relative gap']) if prob[j][2] == 0 else "--", str(sol_2['primal objective']) if prob[j][2] == 0 else "--", \
                                                 str(sol_2['dual objective']) if prob[j][2] == 0 else "--", str(sol_2['primal infeasibility']) if prob[j][2] == 0 else "--", str(sol_2['dual infeasibility']) if prob[j][2] == 0 else "--", \
                                                 str(sol_2['primal slack']) if prob[j][2] == 0 else "--", str(sol_2['dual slack']) if prob[j][2] == 0 else "--"],
                        '1-norm': [str(sol_4['gap']) if prob[j][3] == 0 else "--", str(sol_4['relative gap']) if prob[j][3] == 0 else "--", str(sol_4['primal objective']) if prob[j][3] == 0 else "--", \
                                               str(sol_4['dual objective']) if prob[j][3] == 0 else "--", str(sol_4['primal infeasibility']) if prob[j][3] == 0 else "--", str(sol_4['dual infeasibility']) if prob[j][3] == 0 else "--", \
                                               str(sol_4['primal slack']) if prob[j][3] == 0 else "--", str(sol_4['dual slack']) if prob[j][3] == 0 else "--"]}
                udaje2 = {'': ['Hard', 'Tight hard', '2-norm^2', '1-norm'],
                        'Status': [sol_3['status'] if prob[j][0] == 0 else "--", sol['status'] if prob[j][1] == 0 else "--", sol_2['status'] if prob[j][2] == 0 else "--", sol_4['status'] if prob[j][3] == 0 else "--"]}
                with open("solution_data_" + fileURT.split('.')[0] + model_type[j] + "_M.txt", h_var) as f:
                    if dtX != 1: f.write('\n\n')
                    f.write('\n' + a[4] + str(config_values.columns[dtX].strftime("%Y-%m-%d %H:%M:%S")) + a[4])
                    f.write('\n=============\nSolution data\n=============\n\n')
                    f.write(tabulate(pd.DataFrame(data = udaje2), headers = 'keys', showindex = False))
                    f.write('\n\n')
                    f.write(tabulate(pd.DataFrame(data = udaje), headers = 'keys', showindex = False))
                if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                    print()
                    print("\nSolution data\n=============\n")
                    print(tabulate(udaje2, headers = 'keys'))
                    print()
                    print(tabulate(udaje, headers = 'keys'))
                    print()

                if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                    # Calculation time
                    cas = {'HARD': [elapsed_BEZ_SOFTU],
                        'TIGHT HARD': [elapsed_PUK_POVODNE],
                        '2-norm^2': [elapsed_2_NORMA],
                        '1-norm': [elapsed_1_NORMA]};
                    print("\nElapsed Time\n============\n")
                    print(tabulate(cas, headers = 'keys', floatfmt = ".8f"))
                    print()

            # Check solution status
            statuses = {'HARD BOUNDS real': sol3, 'TIGHT HARD BOUNDS real': sol, '2-norm^2 BOUNDS real': sol2, '1-norm BOUNDS real': sol4, 
                        'HARD BOUNDS ideal': sol_ideal3, 'TIGHT HARD BOUNDS ideal': sol_ideal, '2-norm^2 BOUNDS ideal': sol_ideal2, '1-norm BOUNDS ideal': sol_ideal4}
            descs = ['_tight_hard', '_2_norm', '_hard', '_1_norm', '_tight_hard', '_2_norm', '_hard', '_1_norm']
            for desc, stat in enumerate(statuses.keys()):
                if prob[int(desc/4)][desc%4] == 0:
                    if statuses[stat]['status'] != 'optimal' and (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Solve is True:
                        if 'BFCCU2C1_data' in sys.argv[-2]:
                            print("\n!!!! Solver did not find optimal solution for " + stat + " model !!!!\n\n")
                        else:
                            print('\n' + stylize("Solver did not find optimal solution for " + stat + " model", fg('dark_orange') + attr('bold')) + '\n\n')
                        table_name = 'Recommendations_' + fileURT.split('.')[0] + descs[desc] + '.txt'
                        with open(table_name, h_var) as file:
                            message = '\nThere is no optimal solution. Gap: ' + str(statuses[stat]['gap']) + '\n\n'
                            file.write(message)
                            file.close()
                elif args_data.Solve is True:
                    if 'BFCCU2C1_data' in sys.argv[-2]:
                        print("\n!!!! Solver converged to INF for " + stat + " model !!!!\n\n")
                    else:
                        print('\n' + stylize("Solver converged to INF for " + stat + " model", fg('dark_orange') + attr('bold')) + '\n\n')
                    table_name = 'Recommendations_' + fileURT.split('.')[0] + descs[desc] + '.txt'
                    with open(table_name, h_var) as file:
                        message = '\nThere is no optimal solution. Gap: INF\n\n'
                        file.write(message)
                        file.close()

            # Solutions dicitonaries
            solution[dtX] = {'real': sol, 'ideal': sol_ideal}
            solution2[dtX] = {'real': sol2, 'ideal': sol_ideal2}
            solution3[dtX] = {'real': sol3, 'ideal': sol_ideal3}
            solution4[dtX] = {'real': sol4, 'ideal': sol_ideal4}

            # Insert steady state data
            data_CV[dtX], data_MV[dtX] = insert_solution_into_data(solution[dtX], solution2[dtX], solution3[dtX], solution4[dtX], data_CV[dtX], data_MV[dtX], prob, 2)

        # Calculations solved
        if args_data.Solve is True:
            if 'BFCCU2C1_data' in sys.argv[-2]:
                print('\n!!!! Calculations of optimal point solved !!!!\n\n')
            else:
                print('\n' + stylize('Calculations of optimal point solved', fg('green') + attr('bold')) + '\n\n')
        
        # Save data to external file
        if os.path.isfile("data_" + fileURT + "_" + fileExcel + "_M.txt") is True:
            t1, t2 = load_from_file("data_" + fileURT + "_" + fileExcel + "_M.txt", 0)
            (_, _, _, _, _, _, _, _, mtime_URT, _) = os.stat(fileURT)
            (_, _, _, _, _, _, _, _, mtime_XLSX, _) = os.stat(fileExcel)
            if t1 < int(mtime_URT) or t2 < int(mtime_XLSX):
                save_to_file(solution, solution2, solution3, solution4, M, M2, M3, M4, data_CV, data_MV, data, fileExcel, fileURT)
                if args_data.Solve is True:
                    if 'BFCCU2C1_data' in sys.argv[-2]:
                        print('\n!!!! Data was saved to external file !!!!\n\n')
                    else:
                        print('\n' + stylize('Data was saved to external file', fg('blue') + attr('bold')) + '\n\n')
        else:
            save_to_file(solution, solution2, solution3, solution4, M, M2, M3, M4, data_CV, data_MV, data, fileExcel, fileURT)
            if args_data.Solve is True:
                if 'BFCCU2C1_data' in sys.argv[-2]:
                    print('\n!!!! Data was saved to external file !!!!\n\n')
                else:
                    print('\n' + stylize('Data was saved to external file', fg('blue') + attr('bold')) + '\n\n')

        if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Solve is True:
            os.system('start "" /max ' + 'limits_' + fileURT.split('.')[0] + '_real_M.txt')
            os.system('start "" /max ' + 'limits_' + fileURT.split('.')[0] + '_ideal_M.txt')
            os.system('start "" /max ' + 'objective_' + fileURT.split('.')[0] + '_real_M.txt')
            os.system('start "" /max ' + 'objective_' + fileURT.split('.')[0] + '_ideal_M.txt')
            os.system('start "" /max ' + 'epsilons_2_norm_' + fileURT.split('.')[0] + '_real_M.txt')
            os.system('start "" /max ' + 'epsilons_2_norm_' + fileURT.split('.')[0] + '_ideal_M.txt')
            os.system('start "" /max ' + 'epsilons_1_norm_' + fileURT.split('.')[0] + '_real_M.txt')
            os.system('start "" /max ' + 'epsilons_1_norm_' + fileURT.split('.')[0] + '_ideal_M.txt')
            os.system('start "" /max ' + 'solution_data_' + fileURT.split('.')[0] + '_real_M.txt')
            os.system('start "" /max ' + 'solution_data_' + fileURT.split('.')[0] + '_ideal_M.txt')

    except Exception as e:
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)
        if 'BFCCU2C1_data' in sys.argv[-2]:
            print('\n!!!! Load correct files and solve problem !!!!\n\n')
        else:
            print('\n' + stylize('Load correct files and solve problem', fg('red') + attr('bold')) + '\n\n')

def save_to_file(solution, solution2, solution3, solution4, M, M2, M3, M4, data_CV, data_MV, data, fileExcel, fileURT):
    try:
        if "SteadyStateBFCCU2C1.xlsx" in fileExcel:
            solution, solution2, solution3, solution4, M, M2, M3, M4, data_CV, data_MV, data = \
                {1: solution}, {1: solution2}, {1: solution3}, {1: solution4}, {1: M}, {1: M2}, {1: M3}, {1: M4}, \
                                {1: data_CV}, {1: data_MV}, {1: data}
        solutionN, solution2N, solution3N, solution4N = solution, solution2, solution3, solution4
        Mn, M2n, M3n, M4n = M, M2, M3, M4
        sols = [solutionN, solution2N, solution3N, solution4N]
        Ms = [Mn, M2n, M3n, M4n]

        for i in range(4):
            for j in range(1, len(data.items()) + 1):
                for k in ['real', 'ideal']:
                    if sols[i][j][k] != 0:
                        sols[i][j][k]['s'] = np.array(sols[i][j][k]['s']).tolist()
                        sols[i][j][k]['x'] = np.array(sols[i][j][k]['x']).tolist()
                        sols[i][j][k]['y'] = np.array(sols[i][j][k]['y']).tolist()
                        sols[i][j][k]['z'] = np.array(sols[i][j][k]['z']).tolist()

        for i in range(4):
            for j in range(1, len(data.items()) + 1):
                Ms[i][j]['A'] = np.array(Ms[i][j]['A']).tolist()
                Ms[i][j]['G'] = np.array(Ms[i][j]['G']).tolist()
                Ms[i][j]['P'] = np.array(Ms[i][j]['P']).tolist()
                Ms[i][j]['b'] = np.array(Ms[i][j]['b']).tolist()
                Ms[i][j]['h'] = np.array(Ms[i][j]['h']).tolist()
                Ms[i][j]['h_ideal'] = np.array(Ms[i][j]['h_ideal']).tolist()
                Ms[i][j]['q'] = np.array(Ms[i][j]['q']).tolist()
                Ms[i][j]['r'] = np.array(Ms[i][j]['r']).tolist()

        (_, _, _, _, _, _, _, _, mtime_URT, _) = os.stat(fileURT)
        (_, _, _, _, _, _, _, _, mtime_XLSX, _) = os.stat(fileExcel)
        dataN = {"fileURT_time": int(mtime_URT), \
                "fileExcel_time": int(mtime_XLSX), \
                "solution": solutionN, "solution2": solution2N, "solution3": solution3N, "solution4": solution4N, \
                "M": Mn, "M2": M2n, "M3": M3n, "M4": M4n, \
                "data_CV": data_CV, "data_MV": data_MV, "data": data}

        if "SteadyStateBFCCU2C1.xlsx" in fileExcel:
            with open("data_" + fileURT + "_" + fileExcel + "_S.txt", "w") as d:
                json.dump(dataN, d)
        else:
            with open("data_" + fileURT + "_" + fileExcel + "_M.txt", "w") as d:
                json.dump(dataN, d)
    except Exception as e:
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)


def load_from_file(fileN, idx):
    try:
        with open(fileN) as f:
            data = json.loads(f.read())
        
        fileURT_time, fileExcel_time, solution, solution2, solution3, solution4, M, M2, M3, M4, data_CV, data_MV, data = \
            itemgetter("fileURT_time", "fileExcel_time", "solution", "solution2", "solution3", "solution4", "M", "M2", "M3", "M4", "data_CV", "data_MV", "data")(data)

        if idx == 0:
            return int(fileURT_time), int(fileExcel_time)
        elif idx == 1:
            sols = [solution, solution2, solution3, solution4]
            Ms = [M, M2, M3, M4]
            datas = [data_CV, data_MV]

            for i in range(4):
                sols[i] = {int(u):v for u,v in sols[i].items()}
                for j in range(1, len(data.items()) + 1):
                    for k in ['real', 'ideal']:
                        if sols[i][j][k] != 0:
                            sols[i][j][k]['s'] = matrix(np.array(sols[i][j][k]['s']))
                            sols[i][j][k]['x'] = matrix(np.array(sols[i][j][k]['x']))
                            sols[i][j][k]['y'] = matrix(np.array(sols[i][j][k]['y']))
                            sols[i][j][k]['z'] = matrix(np.array(sols[i][j][k]['z']))
            solution, solution2, solution3, solution4 = sols

            for i in range(4):
                Ms[i] = {int(u):v for u,v in Ms[i].items()}
                for j in range(1, len(data.items()) + 1):
                    Ms[i][j]['A'] = matrix(np.array(Ms[i][j]['A']))
                    Ms[i][j]['G'] = matrix(np.array(Ms[i][j]['G']))
                    Ms[i][j]['P'] = matrix(np.array(Ms[i][j]['P']))
                    Ms[i][j]['b'] = matrix(np.array(Ms[i][j]['b']))
                    Ms[i][j]['h'] = matrix(np.array(Ms[i][j]['h']))
                    Ms[i][j]['h_ideal'] = matrix(np.array(Ms[i][j]['h_ideal']))
                    Ms[i][j]['q'] = matrix(np.array(Ms[i][j]['q']))
                    Ms[i][j]['r'] = matrix(np.array(Ms[i][j]['r']))
            M, M2, M3, M4 = Ms

            for i in range(2):
                datas[i] = {int(u):v for u,v in datas[i].items()}
                for j in range(1, len(data.items()) + 1):
                    datas[i][j] = {int(w):x for w,x in datas[i][j].items()}
            data_CV, data_MV = datas

            data = {int(u):v for u,v in data.items()}
            for i in range(1, len(data.items()) + 1):
                data[i] = {int(u):v for u,v in data[i].items()}
                for j in range(len(data[i])):
                    data[i][j] = {int(w):x for w,x in data[i][j].items()}

            if "SteadyStateBFCCU2C1.xlsx" in fileN:
                return solution[1], solution2[1], solution3[1], solution4[1], M[1], M2[1], M3[1], M4[1], data_CV[1], data_MV[1], data[1]
            else:
                return solution, solution2, solution3, solution4, M, M2, M3, M4, data_CV, data_MV, data
    except Exception as e:
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)

# Calculate loss
def cost_loss(solution, solution2, solution3, solution4, M, M2, M3, M4, total, index, data_CV, data_MV, case = 0):
    try:
        solvers.options['show_progress'] = False
        Tlen = len(data_CV) + len(data_MV)
        solvers.options['verbose'] = False

        # 1-norm bounds
        if case == 0:
            # Modify the right side
            _, _, _, h_moded, _, _ = create_matrix_h(data_CV, data_MV, ideal = 0)

            # Modify the right-hand side
            if total == 0:
                h_moded[index] = M4['h_ideal'][index]
            else:
                h_moded = M4['h_ideal']

            # Solve optimization problem
            sol = 0.0
            sol_loss = solvers.qp(2*M4['P'], M4['q'], M4['G'], h_moded, M4['A'], M4['b'], solver = "mosek")

            # Compare objective function change
            J1 = solution4['real']['x'][0:Tlen, :].trans() * M4['P'][0:Tlen, 0:Tlen] * solution4['real']['x'][0:Tlen, :] + M4['q'][0:Tlen, :].trans() * solution4['real']['x'][0:Tlen, :]
            J2 = sol_loss['x'][0:Tlen, :].trans() * M4['P'][0:Tlen, 0:Tlen] * sol_loss['x'][0:Tlen, :] + M4['q'][0:Tlen, :].trans() * sol_loss['x'][0:Tlen, :]
            loss = float(np.array(J1 - J2))
            Rloss = float(np.array(J1 - J2)/abs(np.array(J1))*100.0)

         # Hard bounds
        elif case == 2:

            # Create matrix h
            h_moded, _, _, h_moded1, _, _ = create_matrix_h(data_CV, data_MV, ideal = 0)

            # Move the border/s by the amount of the penalty
            h_moded[0:len(data_CV), :] = h_moded[0:len(data_CV), :] + matrix(np.array([i if i>1e-2 else 0.0 for i in solution4['real']['x'][Tlen:Tlen + len(data_CV), :]]))
            h_moded[len(data_CV):2*len(data_CV), :] = h_moded[len(data_CV):2*len(data_CV), :] - matrix(np.array([i if i>1e-2 else 0.0 for i in solution4['real']['x'][Tlen + len(data_CV):Tlen + 2*len(data_CV), :]]))
            h_moded[2*len(data_CV):Tlen + len(data_CV), :] = h_moded[2*len(data_CV):Tlen + len(data_CV), :] + matrix(np.array([i if i>1e-2 else 0.0 for i in solution4['real']['x'][Tlen + 2*len(data_CV):2*Tlen + len(data_CV), :]]))
            h_moded[Tlen + len(data_CV):2*Tlen, :] = h_moded[Tlen + len(data_CV):2*Tlen, :] - matrix(np.array([i if i>1e-2 else 0.0 for i in solution4['real']['x'][2*Tlen + len(data_CV):3*Tlen, :]]))

            # Solve optimization problem before shift - Tight hard bounds
            sol = solvers.qp(2*M['P'], M['q'], M['G'], h_moded, M['A'], M['b'], solver = "mosek")

            # Move the border/s back by the amount of the penalty
            h_moded[0:len(data_CV), :] = h_moded[0:len(data_CV), :] - matrix(np.array([i if i>1e-2 else 0.0 for i in solution4['real']['x'][Tlen:Tlen + len(data_CV), :]]))
            h_moded[len(data_CV):2*len(data_CV), :] = h_moded[len(data_CV):2*len(data_CV), :] + matrix(np.array([i if i>1e-2 else 0.0 for i in solution4['real']['x'][Tlen + len(data_CV):Tlen + 2*len(data_CV), :]]))
            h_moded[2*len(data_CV):Tlen + len(data_CV), :] = h_moded[2*len(data_CV):Tlen + len(data_CV), :] - matrix(np.array([i if i>1e-2 else 0.0 for i in solution4['real']['x'][Tlen + 2*len(data_CV):2*Tlen + len(data_CV), :]]))
            h_moded[Tlen + len(data_CV):2*Tlen, :] = h_moded[Tlen + len(data_CV):2*Tlen, :] + matrix(np.array([i if i>1e-2 else 0.0 for i in solution4['real']['x'][2*Tlen + len(data_CV):3*Tlen, :]]))

            # Modify the right-hand side
            if total == 0:
                h_moded[index] = M['h_ideal'][index]
                h_moded1[index] = M4['h_ideal'][index]
            else:
                h_moded = M['h_ideal']
                h_moded1 = M4['h_ideal']

            # Check for feasible solution - 1-norm soft bounds
            sol_1norm = solvers.qp(2*M4['P'], M4['q'], M4['G'], h_moded1, M4['A'], M4['b'], solver = "mosek")

            # Move the border/s by the amount of the penalty
            h_moded[0:len(data_CV), :] = h_moded[0:len(data_CV), :] + matrix(np.array([i if i>1e-2 else 0.0 for i in sol_1norm['x'][Tlen:Tlen + len(data_CV), :]]))
            h_moded[len(data_CV):2*len(data_CV), :] = h_moded[len(data_CV):2*len(data_CV), :] - matrix(np.array([i if i>1e-2 else 0.0 for i in sol_1norm['x'][Tlen + len(data_CV):Tlen + 2*len(data_CV), :]]))
            h_moded[2*len(data_CV):Tlen + len(data_CV), :] = h_moded[2*len(data_CV):Tlen + len(data_CV), :] + matrix(np.array([i if i>1e-2 else 0.0 for i in sol_1norm['x'][Tlen + 2*len(data_CV):2*Tlen + len(data_CV), :]]))
            h_moded[Tlen + len(data_CV):2*Tlen, :] = h_moded[Tlen + len(data_CV):2*Tlen, :] - matrix(np.array([i if i>1e-2 else 0.0 for i in sol_1norm['x'][2*Tlen + len(data_CV):3*Tlen, :]]))

            # Solve optimization problem after shift - Tight hard bounds
            sol_loss = solvers.qp(2*M3['P'], M3['q'], M3['G'], h_moded, M3['A'], M3['b'], solver = "mosek")

            # Compare objective function change
            J1 = sol['x'][0:Tlen, :].trans() * M['P'][0:Tlen, 0:Tlen] * sol['x'][0:Tlen, :] + M['q'][0:Tlen, :].trans() * sol['x'][0:Tlen, :]
            J2 = sol_loss['x'][0:Tlen, :].trans() * M['P'][0:Tlen, 0:Tlen] * sol_loss['x'][0:Tlen, :] + M['q'][0:Tlen, :].trans() * sol_loss['x'][0:Tlen, :]
            loss = float(np.array(J1 - J2))
            Rloss = float(np.array(J1 - J2)/abs(np.array(J1))*100.0)

    except Exception as e:
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)

    return loss, Rloss, J1, J2, sol, sol_loss


# Check the effect of shifting limits - combination of losses
def check_limit_comb(pert, pert1, d_CV, d_MV, solution, solution4, M, M4):
    try:
        pert["ID"], pert1["ID"] = pert.index, pert1.index
        perts = [np.array(pert[pert['Recommendation'].str.contains("Move")]['ID']).tolist(), \
                np.array(pert1[pert1['Recommendation'].str.contains("Move")]['ID']).tolist()]
        sols, Ms, loss, Rloss, old_J, new_J = [solution, solution4], [M, M4], [], [], [], []
        solvers.options['show_progress'] = False
        Tlen = len(d_CV) + len(d_MV)

        # Iterate over pertubation models
        for j, i in enumerate(perts):

            if not i:
                loss.append(0.0), Rloss.append(0.0), old_J.append(0.0), new_J.append(0.0)
                continue

            # Modify the right side
            h_moded_, _, _, h_moded1, _, _ = create_matrix_h(d_CV, d_MV, ideal = 0)
            if j == 0:
                # Move the border/s by the amount of the penalty
                h_moded_[0:len(d_CV), :] = h_moded_[0:len(d_CV), :] + matrix(np.array([i if i>1e-2 else 0.0 for i in solution4['real']['x'][Tlen:Tlen + len(d_CV), :]]))
                h_moded_[len(d_CV):2*len(d_CV), :] = h_moded_[len(d_CV):2*len(d_CV), :] - matrix(np.array([i if i>1e-2 else 0.0 for i in solution4['real']['x'][Tlen + len(d_CV):Tlen + 2*len(d_CV), :]]))
                h_moded_[2*len(d_CV):Tlen + len(d_CV), :] = h_moded_[2*len(d_CV):Tlen + len(d_CV), :] + matrix(np.array([i if i>1e-2 else 0.0 for i in solution4['real']['x'][Tlen + 2*len(d_CV):2*Tlen + len(d_CV), :]]))
                h_moded_[Tlen + len(d_CV):2*Tlen, :] = h_moded_[Tlen + len(d_CV):2*Tlen, :] - matrix(np.array([i if i>1e-2 else 0.0 for i in solution4['real']['x'][2*Tlen + len(d_CV):3*Tlen, :]]))

                # Solve optimization problem before shift - Tight hard bounds
                sols[j]['real'] = solvers.qp(2*M['P'], M['q'], M['G'], h_moded_, M['A'], M['b'], solver = "mosek")

                # Move the border/s back by the amount of the penalty
                h_moded_[0:len(d_CV), :] = h_moded_[0:len(d_CV), :] - matrix(np.array([i if i>1e-2 else 0.0 for i in solution4['real']['x'][Tlen:Tlen + len(d_CV), :]]))
                h_moded_[len(d_CV):2*len(d_CV), :] = h_moded_[len(d_CV):2*len(d_CV), :] + matrix(np.array([i if i>1e-2 else 0.0 for i in solution4['real']['x'][Tlen + len(d_CV):Tlen + 2*len(d_CV), :]]))
                h_moded_[2*len(d_CV):Tlen + len(d_CV), :] = h_moded_[2*len(d_CV):Tlen + len(d_CV), :] - matrix(np.array([i if i>1e-2 else 0.0 for i in solution4['real']['x'][Tlen + 2*len(d_CV):2*Tlen + len(d_CV), :]]))
                h_moded_[Tlen + len(d_CV):2*Tlen, :] = h_moded_[Tlen + len(d_CV):2*Tlen, :] + matrix(np.array([i if i>1e-2 else 0.0 for i in solution4['real']['x'][2*Tlen + len(d_CV):3*Tlen, :]]))

                h_moded_[i] = Ms[j]['h_ideal'][i]
                h_moded1[i] = M4['h_ideal'][i]

                # Check for feasible solution - 1-norm soft bounds
                sol_1norm = solvers.qp(2*M4['P'], M4['q'], M4['G'], h_moded1, M4['A'], M4['b'], solver = "mosek")

                # Move the border/s by the amount of the penalty
                h_moded_[0:len(d_CV), :] = h_moded_[0:len(d_CV), :] + matrix(np.array([i if i>1e-2 else 0.0 for i in sol_1norm['x'][Tlen:Tlen + len(d_CV), :]]))
                h_moded_[len(d_CV):2*len(d_CV), :] = h_moded_[len(d_CV):2*len(d_CV), :] - matrix(np.array([i if i>1e-2 else 0.0 for i in sol_1norm['x'][Tlen + len(d_CV):Tlen + 2*len(d_CV), :]]))
                h_moded_[2*len(d_CV):Tlen + len(d_CV), :] = h_moded_[2*len(d_CV):Tlen + len(d_CV), :] + matrix(np.array([i if i>1e-2 else 0.0 for i in sol_1norm['x'][Tlen + 2*len(d_CV):2*Tlen + len(d_CV), :]]))
                h_moded_[Tlen + len(d_CV):2*Tlen, :] = h_moded_[Tlen + len(d_CV):2*Tlen, :] - matrix(np.array([i if i>1e-2 else 0.0 for i in sol_1norm['x'][2*Tlen + len(d_CV):3*Tlen, :]]))
                
                h_moded = h_moded_
            elif j == 1:
                h_moded1[i] = Ms[j]['h_ideal'][i]
                h_moded = h_moded1

            # Solve optimization problem
            sol_loss = solvers.qp(2*Ms[j]['P'], Ms[j]['q'], Ms[j]['G'], h_moded, Ms[j]['A'], Ms[j]['b'], solver = "mosek")

            # Compare objective function change
            J1 = sols[j]['real']['x'][0:Tlen, :].trans() * Ms[j]['P'][0:Tlen, 0:Tlen] * sols[j]['real']['x'][0:Tlen, :] + Ms[j]['q'][0:Tlen, :].trans() * sols[j]['real']['x'][0:Tlen, :]
            J2 = sol_loss['x'][0:Tlen, :].trans() * Ms[j]['P'][0:Tlen, 0:Tlen] * sol_loss['x'][0:Tlen, :] + Ms[j]['q'][0:Tlen, :].trans() * sol_loss['x'][0:Tlen, :]
            loss.append(float(abs(np.array(J1 - J2))))
            Rloss.append(float(np.array(J1 - J2)/abs(np.array(J1))*100.0))
            old_J.append(float(np.array(J1))), new_J.append(float(np.array(J2)))

        # Append to the tables
        pert = pd.concat([pd.DataFrame([{"Priority": 0, "Name": "Total loss", "Variable": " ", "Limit": " ", "Relative loss": Rloss[0], "Recommendation": " ", \
                                   "|": "|", "A/i ABS Real": " ", "A/i REL Real": " ", "A/i ABS Ideal": " ", "A/i REL Ideal": " ", "Absolute loss": loss[0], "ID": -1, \
                                    "Old J": old_J[0], "New J": new_J[0]}]), pert], ignore_index = True).drop(columns = ["ID"])
        pert1 = pd.concat([pd.DataFrame([{"Priority": 0, "Name": "Total loss", "Variable": " ", "Limit": " ", "Relative loss": Rloss[1], "Recommendation": " ", \
                                   "|": "|", "A/i ABS Real": " ", "A/i REL Real": " ", "A/i ABS Ideal": " ", "A/i REL Ideal": " ", "Absolute loss": loss[1], "ID": -1, \
                                    "Old J": old_J[1], "New J": new_J[1]}]), pert1], ignore_index = True).drop(columns = ["ID"])

        return pert, pert1

    except Exception as e:
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)


# Check the effect of shifting limits
def check_limit(solution, solution2, solution3, solution4, M, M2, M3, M4, data_CV, data_MV, args_data):

    try:
        CV_len, MV_len, CMV_len = len(data_CV), len(data_MV), len(data_CV) + len(data_MV)

        name, variable, limit, loss, Rloss, recommendation, abs_real, abs_ideal, rel_real, rel_ideal, old_J, new_J = [], [], [], [], [], [], [], [], [], [], [], []
        # Hard bounds, CV, HL
        for i1 in range(CV_len):
            name.append(data_CV[list(data_CV.keys())[i1]]['CV Name'])
            variable.append("CV " + str(list(data_CV.keys())[i1] + 1))
            limit.append("HL")
            help_loss, help_Rloss, help_old_J, help_new_J, _, _ = cost_loss(solution, solution2, solution3, solution4, M, M2, M3, M4, 0, i1, data_CV, data_MV, case = 2)
            loss.append(help_loss), Rloss.append(help_Rloss), old_J.append(float(np.array(help_old_J))), new_J.append(float(np.array(help_new_J)))
            if help_loss > args_data.Atol:
                recommendation.append("Move HL to IDHL.")
            else:
                recommendation.append("Do not move HL.")
            real_s = data_CV[list(data_CV.keys())[i1]]['High Limit'] - float(data_CV[list(data_CV.keys())[i1]]['Delta Soft High Limit']) + solution4['real']['x'][CMV_len + i1]
            ideal_s = data_CV[list(data_CV.keys())[i1]]['Ideal High'] - float(data_CV[list(data_CV.keys())[i1]]['Delta Soft High Limit']) + solution4['ideal']['x'][CMV_len + i1]
            if abs(solution4['real']['x'][i1] - (real_s)) <= args_data.Acomp:
                abs_real.append('ACTIVE')
            else:
                abs_real.append('inactive')
            if abs(solution4['ideal']['x'][i1] - (ideal_s)) <= args_data.Acomp:
                abs_ideal.append('ACTIVE')
            else:
                abs_ideal.append('inactive')
            try:
                if abs(solution4['real']['x'][i1] - (real_s))/real_s <= args_data.Rcomp:
                    rel_real.append('ACTIVE')
                else:
                    rel_real.append('inactive')
            except:
                rel_real.append(abs_real[-1])
            try:
                if abs(solution4['ideal']['x'][i1] - (ideal_s))/ideal_s <= args_data.Rcomp:
                    rel_ideal.append('ACTIVE')
                else:
                    rel_ideal.append('inactive')
            except:
                rel_ideal.append(abs_ideal[-1])

        # Hard bounds, CV, LL
        for i1 in range(CV_len):
            i2 = i1 + CV_len
            name.append(data_CV[list(data_CV.keys())[i1]]['CV Name'])
            variable.append("CV " + str(list(data_CV.keys())[i1] + 1))
            limit.append("LL")
            help_loss, help_Rloss, help_old_J, help_new_J, _, _ = cost_loss(solution, solution2, solution3, solution4, M, M2, M3, M4, 0, i2, data_CV, data_MV, case = 2)
            loss.append(help_loss), Rloss.append(help_Rloss), old_J.append(float(np.array(help_old_J))), new_J.append(float(np.array(help_new_J)))
            if help_loss > args_data.Atol:
                recommendation.append("Move LL to IDLL.")
            else:
                recommendation.append("Do not move LL.")
            real_s = data_CV[list(data_CV.keys())[i1]]['Low Limit'] + float(data_CV[list(data_CV.keys())[i1]]['Delta Soft Low Limit'] )- solution4['real']['x'][CMV_len + i2]
            ideal_s = data_CV[list(data_CV.keys())[i1]]['Ideal Low'] + float(data_CV[list(data_CV.keys())[i1]]['Delta Soft Low Limit'] )- solution4['ideal']['x'][CMV_len + i2]
            if abs(solution4['real']['x'][i1] - (real_s)) <= args_data.Acomp:
                abs_real.append('ACTIVE')
            else:
                abs_real.append('inactive')
            if abs(solution4['ideal']['x'][i1] - (ideal_s)) <= args_data.Acomp:
                abs_ideal.append('ACTIVE')
            else:
                abs_ideal.append('inactive')
            try:
                if abs(solution4['real']['x'][i1] - (real_s))/real_s <= args_data.Rcomp:
                    rel_real.append('ACTIVE')
                else:
                    rel_real.append('inactive')
            except:
                rel_real.append(abs_real[-1])
            try:
                if abs(solution4['ideal']['x'][i1] - (ideal_s))/ideal_s <= args_data.Rcomp:
                    rel_ideal.append('ACTIVE')
                else:
                    rel_ideal.append('inactive')
            except:
                rel_ideal.append(abs_ideal[-1])

        # Hard bounds, MV, HL
        for i1 in range(MV_len):
            i2 = i1 + 2*CV_len
            i3 = i1 + CV_len
            name.append(data_MV[list(data_MV.keys())[i1]]['MV Name'])
            variable.append("MV " + str(list(data_MV.keys())[i1] + 1))
            limit.append("HL")
            help_loss, help_Rloss, help_old_J, help_new_J, _, _ = cost_loss(solution, solution2, solution3, solution4, M, M2, M3, M4, 0, i2, data_CV, data_MV, case = 2)
            loss.append(help_loss), Rloss.append(help_Rloss), old_J.append(float(np.array(help_old_J))), new_J.append(float(np.array(help_new_J)))
            if help_loss > args_data.Atol:
                recommendation.append("Move HL to IDHL.")
            else:
                recommendation.append("Do not move HL.")
            real_s = data_MV[list(data_MV.keys())[i1]]['High Limit'] - float(data_MV[list(data_MV.keys())[i1]]['Delta Soft High Limit']) + solution4['real']['x'][CMV_len + i2]
            ideal_s = data_MV[list(data_MV.keys())[i1]]['Ideal High'] - float(data_MV[list(data_MV.keys())[i1]]['Delta Soft High Limit']) + solution4['ideal']['x'][CMV_len + i2]
            if abs(solution4['real']['x'][i3] - (real_s)) <= args_data.Acomp:
                abs_real.append('ACTIVE')
            else:
                abs_real.append('inactive')
            if abs(solution4['ideal']['x'][i3] - (ideal_s)) <= args_data.Acomp:
                abs_ideal.append('ACTIVE')
            else:
                abs_ideal.append('inactive')
            try:
                if abs(solution4['real']['x'][i3] - (real_s))/real_s <= args_data.Rcomp:
                    rel_real.append('ACTIVE')
                else:
                    rel_real.append('inactive')
            except:
                rel_real.append(abs_real[-1])
            try:
                if abs(solution4['ideal']['x'][i3] - (ideal_s))/ideal_s <= args_data.Rcomp:
                    rel_ideal.append('ACTIVE')
                else:
                    rel_ideal.append('inactive')
            except:
                rel_ideal.append(abs_ideal[-1])

        # Hard bounds, MV, LL
        for i1 in range(MV_len):
            i2 = i1 + 2*CV_len + MV_len
            i3 = i1 + CV_len
            name.append(data_MV[list(data_MV.keys())[i1]]['MV Name'])
            variable.append("MV " + str(list(data_MV.keys())[i1] + 1))
            limit.append("LL")
            help_loss, help_Rloss, help_old_J, help_new_J, _, _ = cost_loss(solution, solution2, solution3, solution4, M, M2, M3, M4, 0, i2, data_CV, data_MV, case = 2)
            loss.append(help_loss), Rloss.append(help_Rloss), old_J.append(float(np.array(help_old_J))), new_J.append(float(np.array(help_new_J)))
            if help_loss > args_data.Atol:
                recommendation.append("Move LL to IDLL.")
            else:
                recommendation.append("Do not move LL.")
            real_s = data_MV[list(data_MV.keys())[i1]]['Low Limit'] + float(data_MV[list(data_MV.keys())[i1]]['Delta Soft Low Limit'] )- solution4['real']['x'][CMV_len + i2]
            ideal_s = data_MV[list(data_MV.keys())[i1]]['Ideal Low'] + float(data_MV[list(data_MV.keys())[i1]]['Delta Soft Low Limit'] )- solution4['ideal']['x'][CMV_len + i2]
            if abs(solution4['real']['x'][i3] - (real_s)) <= args_data.Acomp:
                abs_real.append('ACTIVE')
            else:
                abs_real.append('inactive')
            if abs(solution4['ideal']['x'][i3] - (ideal_s)) <= args_data.Acomp:
                abs_ideal.append('ACTIVE')
            else:
                abs_ideal.append('inactive')
            try:
                if abs(solution4['real']['x'][i3] - (real_s))/real_s <= args_data.Rcomp:
                    rel_real.append('ACTIVE')
                else:
                    rel_real.append('inactive')
            except:
                rel_real.append(abs_real[-1])
            try:
                if abs(solution4['ideal']['x'][i3] - (ideal_s))/ideal_s <= args_data.Rcomp:
                    rel_ideal.append('ACTIVE')
                else:
                    rel_ideal.append('inactive')
            except:
                rel_ideal.append(abs_ideal[-1])

        # Table data - Hard bounds
        table = {'Name': name, 'Variable': variable, 'Limit': limit, 'Relative loss': Rloss, 'Recommendation': recommendation, '|': ['|']*len(name), \
                  'A/i ABS Real': abs_real, 'A/i REL Real': rel_real, 'A/i ABS Ideal': abs_ideal, 'A/i REL Ideal': rel_ideal, 'Absolute loss': loss, \
                  'Old J': old_J, 'New J': new_J}
        df = pd.DataFrame(data = table).sort_values(by = ['Relative loss', 'Variable', 'Limit'], ascending = [False, True, True])
        array = np.arange(1, len(name) + 1, dtype = int)
        df.insert(0, "Priority", array)

        name, variable, limit, loss, Rloss, recommendation, abs_real, abs_ideal, rel_real, rel_ideal, old_J, new_J = [], [], [], [], [], [], [], [], [], [], [], []
        # 1-norm, CV, HL
        for i1 in range(CV_len):
            name.append(data_CV[list(data_CV.keys())[i1]]['CV Name'])
            variable.append("CV " + str(list(data_CV.keys())[i1] + 1))
            limit.append("HL")
            help_loss, help_Rloss, help_old_J, help_new_J, _, _ = cost_loss(solution, solution2, solution3, solution4, M, M2, M3, M4, 0, i1, data_CV, data_MV, case = 0)
            loss.append(help_loss), Rloss.append(help_Rloss), old_J.append(float(np.array(help_old_J))), new_J.append(float(np.array(help_new_J)))
            if help_loss > args_data.Atol:
                recommendation.append("Move HL to IDHL.")
            else:
                recommendation.append("Do not move HL.")
            real_s = data_CV[list(data_CV.keys())[i1]]['High Limit'] - float(data_CV[list(data_CV.keys())[i1]]['Delta Soft High Limit']) + solution4['real']['x'][CMV_len + i1]
            ideal_s = data_CV[list(data_CV.keys())[i1]]['Ideal High'] - float(data_CV[list(data_CV.keys())[i1]]['Delta Soft High Limit']) + solution4['ideal']['x'][CMV_len + i1]
            if abs(solution4['real']['x'][i1] - (real_s)) <= args_data.Acomp:
                abs_real.append('ACTIVE')
            else:
                abs_real.append('inactive')
            if abs(solution4['ideal']['x'][i1] - (ideal_s)) <= args_data.Acomp:
                abs_ideal.append('ACTIVE')
            else:
                abs_ideal.append('inactive')
            try:
                if abs(solution4['real']['x'][i1] - (real_s))/real_s <= args_data.Rcomp:
                    rel_real.append('ACTIVE')
                else:
                    rel_real.append('inactive')
            except:
                rel_real.append(abs_real[-1])
            try:
                if abs(solution4['ideal']['x'][i1] - (ideal_s))/ideal_s <= args_data.Rcomp:
                    rel_ideal.append('ACTIVE')
                else:
                    rel_ideal.append('inactive')
            except:
                rel_ideal.append(abs_ideal[-1])

        # 1-norm, CV, LL
        for i1 in range(CV_len):
            i2 = i1 + CV_len
            name.append(data_CV[list(data_CV.keys())[i1]]['CV Name'])
            variable.append("CV " + str(list(data_CV.keys())[i1] + 1))
            limit.append("LL")
            help_loss, help_Rloss, help_old_J, help_new_J, _, _ = cost_loss(solution, solution2, solution3, solution4, M, M2, M3, M4, 0, i2, data_CV, data_MV, case = 0)
            loss.append(help_loss), Rloss.append(help_Rloss), old_J.append(float(np.array(help_old_J))), new_J.append(float(np.array(help_new_J)))
            if help_loss > args_data.Atol:
                recommendation.append("Move LL to IDLL.")
            else:
                recommendation.append("Do not move LL.")
            real_s = data_CV[list(data_CV.keys())[i1]]['Low Limit'] + float(data_CV[list(data_CV.keys())[i1]]['Delta Soft Low Limit'] )- solution4['real']['x'][CMV_len + i2]
            ideal_s = data_CV[list(data_CV.keys())[i1]]['Ideal Low'] + float(data_CV[list(data_CV.keys())[i1]]['Delta Soft Low Limit'] )- solution4['ideal']['x'][CMV_len + i2]
            if abs(solution4['real']['x'][i1] - (real_s)) <= args_data.Acomp:
                abs_real.append('ACTIVE')
            else:
                abs_real.append('inactive')
            if abs(solution4['ideal']['x'][i1] - (ideal_s)) <= args_data.Acomp:
                abs_ideal.append('ACTIVE')
            else:
                abs_ideal.append('inactive')
            try:
                if abs(solution4['real']['x'][i1] - (real_s))/real_s <= args_data.Rcomp:
                    rel_real.append('ACTIVE')
                else:
                    rel_real.append('inactive')
            except:
                rel_real.append(abs_real[-1])
            try:
                if abs(solution4['ideal']['x'][i1] - (ideal_s))/ideal_s <= args_data.Rcomp:
                    rel_ideal.append('ACTIVE')
                else:
                    rel_ideal.append('inactive')
            except:
                rel_ideal.append(abs_ideal[-1])

        # 1-norm, MV, HL
        for i1 in range(MV_len):
            i2 = i1 + 2*CV_len
            i3 = i1 + CV_len
            name.append(data_MV[list(data_MV.keys())[i1]]['MV Name'])
            variable.append("MV " + str(list(data_MV.keys())[i1] + 1))
            limit.append("HL")
            help_loss, help_Rloss, help_old_J, help_new_J, _, _ = cost_loss(solution, solution2, solution3, solution4, M, M2, M3, M4, 0, i2, data_CV, data_MV, case = 0)
            loss.append(help_loss), Rloss.append(help_Rloss), old_J.append(float(np.array(help_old_J))), new_J.append(float(np.array(help_new_J)))
            if help_loss > args_data.Atol:
                recommendation.append("Move HL to IDHL.")
            else:
                recommendation.append("Do not move HL.")
            real_s = data_MV[list(data_MV.keys())[i1]]['High Limit'] - float(data_MV[list(data_MV.keys())[i1]]['Delta Soft High Limit']) + solution4['real']['x'][CMV_len + i2]
            ideal_s = data_MV[list(data_MV.keys())[i1]]['Ideal High'] - float(data_MV[list(data_MV.keys())[i1]]['Delta Soft High Limit']) + solution4['ideal']['x'][CMV_len + i2]
            if abs(solution4['real']['x'][i3] - (real_s)) <= args_data.Acomp:
                abs_real.append('ACTIVE')
            else:
                abs_real.append('inactive')
            if abs(solution4['ideal']['x'][i3] - (ideal_s)) <= args_data.Acomp:
                abs_ideal.append('ACTIVE')
            else:
                abs_ideal.append('inactive')
            try:
                if abs(solution4['real']['x'][i3] - (real_s))/real_s <= args_data.Rcomp:
                    rel_real.append('ACTIVE')
                else:
                    rel_real.append('inactive')
            except:
                rel_real.append(abs_real[-1])
            try:
                if abs(solution4['ideal']['x'][i3] - (ideal_s))/ideal_s <= args_data.Rcomp:
                    rel_ideal.append('ACTIVE')
                else:
                    rel_ideal.append('inactive')
            except:
                rel_ideal.append(abs_ideal[-1])

        # 1-norm, MV, LL
        for i1 in range(MV_len):
            i2 = i1 + 2*CV_len + MV_len
            i3 = i1 + CV_len
            name.append(data_MV[list(data_MV.keys())[i1]]['MV Name'])
            variable.append("MV " + str(list(data_MV.keys())[i1] + 1))
            limit.append("LL")
            help_loss, help_Rloss, help_old_J, help_new_J, _, _ = cost_loss(solution, solution2, solution3, solution4, M, M2, M3, M4, 0, i2, data_CV, data_MV, case = 0)
            loss.append(help_loss), Rloss.append(help_Rloss), old_J.append(float(np.array(help_old_J))), new_J.append(float(np.array(help_new_J)))
            if help_loss > args_data.Atol:
                recommendation.append("Move LL to IDLL.")
            else:
                recommendation.append("Do not move LL.")
            real_s = data_MV[list(data_MV.keys())[i1]]['Low Limit'] + float(data_MV[list(data_MV.keys())[i1]]['Delta Soft Low Limit'] )- solution4['real']['x'][CMV_len + i2]
            ideal_s = data_MV[list(data_MV.keys())[i1]]['Ideal Low'] + float(data_MV[list(data_MV.keys())[i1]]['Delta Soft Low Limit'] )- solution4['ideal']['x'][CMV_len + i2]
            if abs(solution4['real']['x'][i3] - (real_s)) <= args_data.Acomp:
                abs_real.append('ACTIVE')
            else:
                abs_real.append('inactive')
            if abs(solution4['ideal']['x'][i3] - (ideal_s)) <= args_data.Acomp:
                abs_ideal.append('ACTIVE')
            else:
                abs_ideal.append('inactive')
            try:
                if abs(solution4['real']['x'][i3] - (real_s))/real_s <= args_data.Rcomp:
                    rel_real.append('ACTIVE')
                else:
                    rel_real.append('inactive')
            except:
                rel_real.append(abs_real[-1])
            try:
                if abs(solution4['ideal']['x'][i3] - (ideal_s))/ideal_s <= args_data.Rcomp:
                    rel_ideal.append('ACTIVE')
                else:
                    rel_ideal.append('inactive')
            except:
                rel_ideal.append(abs_ideal[-1])

        # Table data - 1-norm
        table1 = {'Name': name, 'Variable': variable, 'Limit': limit, 'Relative loss': Rloss, 'Recommendation': recommendation, '|': ['|']*len(name), \
                  'A/i ABS Real': abs_real, 'A/i REL Real': rel_real, 'A/i ABS Ideal': abs_ideal, 'A/i REL Ideal': rel_ideal, 'Absolute loss': loss, \
                  'Old J': old_J, 'New J': new_J}
        df1 = pd.DataFrame(data = table1).sort_values(by = ['Relative loss', 'Variable', 'Limit'], ascending = [False, True, True])
        array1 = np.arange(1, len(name) + 1, dtype = int)
        df1.insert(0, "Priority", array1)

    except Exception as e:
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)

    return df, df1


# Check, if some soft constraints are violated
def check_violation(solution2, solution4, data_CV, data_MV):
    CV_size, MV_size = int(len(data_CV)), int(len(data_MV))
    var, name, val1, type1, lim1, rec1 = [], [], [], [], [], []

    # Check controlled variables
    for h,i in enumerate(data_CV.keys()):
        skip11, skip12 = 0, 0

        # Variable number & name
        var.append('CV ' + data_CV[i]['CV Index'])
        var.append('CV ' + data_CV[i]['CV Index'])
        name.append(data_CV[i]['CV Name'])
        name.append(data_CV[i]['CV Name'])

        # CV real high limits - 1-norm
        if float(np.array(solution4['real']['x'][h + CV_size + MV_size, :])) > 1.0e-5:
            skip11 = 1
            val1.append(float(np.array(solution4['real']['x'][h + CV_size + MV_size, :])))
            type1.append('REAL')
            lim1.append('HL')
            if data_CV[i]['SS Ideal 1-norm'] >= data_CV[i]['SS Real 1-norm']:
                rec1.append('Loosen constraint by Violation value towards ideal high limit')
            else:
                rec1.append('Loosen constraint by Violation value away from ideal high limit')
        # CV real low limits - 1-norm
        elif float(np.array(solution4['real']['x'][h + 2*CV_size + MV_size, :])) > 1.0e-5:
            skip12 = 1
            val1.append(float(np.array(solution4['real']['x'][h + 2*CV_size + MV_size, :])))
            type1.append('REAL')
            lim1.append('LL')
            if data_CV[i]['SS Ideal 1-norm'] <= data_CV[i]['SS Real 1-norm']:
                rec1.append('Loosen constraint by Violation value towards ideal low limit')
            else:
                rec1.append('Loosen constraint by Violation value away from ideal low limit')
        # CV ideal high limits - 1-norm
        elif float(np.array(solution4['ideal']['x'][h + CV_size + MV_size, :])) > 1.0e-5:
            skip11 = 1
            val1.append(float(np.array(solution4['ideal']['x'][h + CV_size + MV_size, :])))
            type1.append('IDEAL')
            lim1.append('HL')
            if data_CV[i]['SS Real 1-norm'] >= data_CV[i]['SS Ideal 1-norm']:
                rec1.append('Loosen constraint by Violation value towards high limit')
            else:
                rec1.append('Loosen constraint by Violation value away from high limit')
        # CV ideal low limits - 1-norm
        elif float(np.array(solution4['ideal']['x'][h + 2*CV_size + MV_size, :])) > 1.0e-5:
            skip12 = 1
            val1.append(float(np.array(solution4['ideal']['x'][h + 2*CV_size + MV_size, :])))
            type1.append('IDEAL')
            lim1.append('LL')
            if data_CV[i]['SS Real 1-norm'] <= data_CV[i]['SS Ideal 1-norm']:
                rec1.append('Loosen constraint by Violation value towards low limit')
            else:
                rec1.append('Loosen constraint by Violation value away from low limit')

        # Not violated CV constraints - 1-norm
        if skip11 == 0:
            val1.append(0.0)
            type1.append('IDEAL/REAL')
            lim1.append('HL')
            rec1.append('Constraint is not violated')
        if skip12 == 0:
            val1.append(0.0)
            type1.append('IDEAL/REAL')
            lim1.append('LL')
            rec1.append('Constraint is not violated')

    # Check manipulated variables
    for h,i in enumerate(data_MV.keys()):
        skip11, skip12, skip21, skip22 = 0, 0, 0, 0

        var.append('MV ' + data_MV[i]['MV Index'])
        var.append('MV ' + data_MV[i]['MV Index'])
        name.append(data_MV[i]['MV Name'])
        name.append(data_MV[i]['MV Name'])

        # MV real high limits - 1-norm
        if float(np.array(solution4['real']['x'][h + 3*CV_size + MV_size, :])) > 1.0e-5:
            skip11 = 1
            val1.append(float(np.array(solution4['real']['x'][h + 3*CV_size + MV_size, :])))
            type1.append('REAL')
            lim1.append('HL')
            if data_MV[i]['SS Ideal 1-norm'] >= data_MV[i]['SS Real 1-norm']:
                rec1.append('Loosen constraint by Violation value towards ideal high limit')
            else:
                rec1.append('Loosen constraint by Violation value away from ideal high limit')
        # MV real low limits - 1-norm
        elif float(np.array(solution4['real']['x'][h + 3*CV_size + 2*MV_size, :])) > 1.0e-5:
            skip12 = 1
            val1.append(float(np.array(solution4['real']['x'][h + 3*CV_size + 2*MV_size, :])))
            type1.append('REAL')
            lim1.append('LL')
            if data_MV[i]['SS Ideal 1-norm'] <= data_MV[i]['SS Real 1-norm']:
                rec1.append('Loosen constraint by Violation value towards ideal low limit')
            else:
                rec1.append('Loosen constraint by Violation value away from ideal low limit')
        # MV ideal high limits - 1-norm
        elif float(np.array(solution4['ideal']['x'][h + 3*CV_size + MV_size, :])) > 1.0e-5:
            skip11 = 1
            val1.append(float(np.array(solution4['ideal']['x'][h + 3*CV_size + MV_size, :])))
            type1.append('IDEAL')
            lim1.append('HL')
            if data_MV[i]['SS Real 1-norm'] >= data_MV[i]['SS Ideal 1-norm']:
                rec1.append('Loosen constraint by Violation value towards high limit')
            else:
                rec1.append('Loosen constraint by Violation value away from high limit')
        # MV ideal low limits - 1-norm
        elif float(np.array(solution4['ideal']['x'][h + 3*CV_size + 2*MV_size, :])) > 1.0e-5:
            skip12 = 1
            val1.append(float(np.array(solution4['ideal']['x'][h + 3*CV_size + 2*MV_size, :])))
            type1.append('IDEAL')
            lim1.append('LL')
            if data_MV[i]['SS Real 1-norm'] <= data_MV[i]['SS Ideal 1-norm']:
                rec1.append('Loosen constraint by Violation value towards low limit')
            else:
                rec1.append('Loosen constraint by Violation value away from low limit')
        
        # Not violated MV constraints - 1-norm
        if skip11 == 0:
            val1.append(0.0)
            type1.append('IDEAL/REAL')
            lim1.append('HL')
            rec1.append('Constraint is not violated')
        if skip12 == 0:
            val1.append(0.0)
            type1.append('IDEAL/REAL')
            lim1.append('LL')
            rec1.append('Constraint is not violated')

    # Table data - 1-norm
    table1 = {'Name': name, 'Variable': var, 'Limit': lim1, 'Model': type1, 'Violation': val1, 'Recommendation': rec1}
    td1 = pd.DataFrame(data = table1).sort_values(by = ['Violation','Variable'], ascending = [False,True])
    array1 = np.arange(1, len(name) + 1, dtype = int)
    td1.insert(0, "Priority", array1)

    return td1


# Calculate non-local analysis
def perturbations_analysis(P, q, G, A, b, sol, case, num, data_CV, data_MV, args_data, solution4, M4):
    eps, e_load = args_data.Change, 0
    CV_size, MV_size = len(data_CV), len(data_MV)
    Tlen = CV_size + MV_size
    s = 2*Tlen
    if case == 0:
        lagra = np.zeros((s, 4))
    else:
        lagra = np.zeros((s, 5))

    # Original objective function (without penalization of 'e' and 'E' parameters)
    if num == 0:
        h_mode, _, _, _, _, _ = create_matrix_h(data_CV, data_MV, ideal = 0)
        h_mode[0:len(data_CV), :] = h_mode[0:len(data_CV), :] + matrix(np.array([i+0.0 if i>1e-2 else 0.0 for i in solution4['real']['x'][Tlen:Tlen + len(data_CV), :]]))
        h_mode[len(data_CV):2*len(data_CV), :] = h_mode[len(data_CV):2*len(data_CV), :] - matrix(np.array([i+0.0 if i>1e-2 else 0.0 for i in solution4['real']['x'][Tlen + len(data_CV):Tlen + 2*len(data_CV), :]]))
        h_mode[2*len(data_CV):Tlen + len(data_CV), :] = h_mode[2*len(data_CV):Tlen + len(data_CV), :] + matrix(np.array([i+0.0 if i>1e-2 else 0.0 for i in solution4['real']['x'][Tlen + 2*len(data_CV):2*Tlen + len(data_CV), :]]))
        h_mode[Tlen + len(data_CV):2*Tlen, :] = h_mode[Tlen + len(data_CV):2*Tlen, :] - matrix(np.array([i+0.0 if i>1e-2 else 0.0 for i in solution4['real']['x'][2*Tlen + len(data_CV):3*Tlen, :]]))
        h_copy = h_mode
        sol = solvers.qp(P, q, G, h_mode, A, b)
    else:
        _, _, _, h_mode, _, _ = create_matrix_h(data_CV, data_MV, ideal = 0)
        h_copy = h_mode

    lm = sol['z'][:s]

    if sol['status'] != 'optimal':
        # Error message

        if num == 1:
            if 'BFCCU2C1_data' in sys.argv[-2]:
                print('\n!!!! Solver status is not optimal - 1-norm soft bounds perturbation analysis might by inaccurate !!!!')
            else:
                print('\n' + stylize('Solver status is not optimal - 1-norm soft bounds perturbation analysis might by inaccurate', fg('dark_orange') + attr('bold')))
        else:
            if 'BFCCU2C1_data' in sys.argv[-2]:
                print('\n!!!! Solver status is not optimal - hard bounds perturbation analysis might by inaccurate !!!!')
            else:
                print('\n' + stylize('Solver status is not optimal - hard bounds perturbation analysis might by inaccurate', fg('dark_orange') + attr('bold')))

    PO1 = sol['primal objective']

    # Round Lagrange multipliers
    lagra[:,1] = np.array([round(i, 4) if i >= 0.1 else round(i, 0) for i in lm])

    # Non-local analysis
    for i in range(s):

        # Load ideal 'h' matrix
        if case == 1:
            if num == 1:
                _, _, _, h_mode_ideal, _, _ = create_matrix_h(data_CV, data_MV, ideal = 1)
            else:
                h_mode_ideal, _, _, _, _, _ = create_matrix_h(data_CV, data_MV, ideal = 1)

            # Update the change of the right side of the inequality constraints
            eps = abs(h_mode_ideal[i] - h_mode[i])
            lagra[i][4] = eps
        else:
            eps = args_data.Change

        h_copy[i] += eps

        # New objective function (without penalization of 'e' and 'E' parameters)
        try:
            sol2 = solvers.qp(P, q, G, h_copy, A, b)

            # New objective function (without penalization of 'e' and 'E' parameters)
            PO2 = sol2['primal objective']

        except:
            PO2 = PO1

        h_copy[i] -= eps

        # New Lagrange multipliers
        try:
            lagra[i][0] = abs(np.array(PO2 - PO1))/eps + 0.0 if abs(np.array(PO2 - PO1))/eps + 0.0 >= 0.1 else 0.0
        except:
            lagra[i][0] = lagra[i][1]

        # Perturbations                
        try:
            lagra[i][2] = lagra[i][0] - lagra[i][1]
            lagra[i][3] = lagra[i][2]/abs(lagra[i][1])*100
        except RuntimeWarning:
            lagra[i][3] = 0.0

        # Update the change of the right side of the inequality constraints
        eps = args_data.Change

    # Show error
    if e_load > 0:
        if num == 1:
            if 'BFCCU2C1_data' in sys.argv[-2]:
                print('\n!!!! Some low constraint/s (' + str(e_load) + ' in total) has higher value than high constraint - 1-norm soft bounds !!!!\n\n')
            else:
                print('\n' + stylize('Some low constraint/s (' + str(e_load) + ' in total) has higher value than high constraint - 1-norm soft bounds', fg('orange') + attr('bold')) + '\n\n')
        if num == 2:
            if 'BFCCU2C1_data' in sys.argv[-2]:
                print('\n!!!! Some low constraint/s (' + str(e_load) + ' in total) has higher value than high constraint - 2-norm^2 soft bounds !!!!\n\n')
            else:
                print('\n' + stylize('Some low constraint/s (' + str(e_load) + ' in total) has higher value than high constraint - 2-norm^2 soft bounds', fg('orange') + attr('bold')) + '\n\n')
        else:
            if 'BFCCU2C1_data' in sys.argv[-2]:
                print('\n!!!! Some low constraint/s (' + str(e_load) + ' in total) has higher value than high constraint - hard bounds !!!!\n\n')
            else:
                print('\n' + stylize('Some low constraint/s (' + str(e_load) + ' in total) has higher value than high constraint - hard bounds', fg('orange') + attr('bold')) + '\n\n')

    return np.around(lagra, 6)


# Non-local analysis - table data
def table_perturbations(x, data_CV, data_MV, model):
    table = pd.DataFrame(x)
    if model == 0:
        table.columns = ['Perturbation', 'Lagrange', 'Absolute error', 'Relative error']
    else:
        table.columns = ['Perturbation', 'Lagrange', 'Absolute error', 'Relative error', 'Shift size']

    CV_name, CV_var, MV_name, MV_var = [], [], [], []

    # CV description
    for i in data_CV.keys():
        CV_name.append(data_CV[i]['CV Name'])
        CV_var.append('CV ' + str(data_CV[i]['CV Index']))

    # MV description
    for i in data_MV.keys():
        MV_name.append(data_MV[i]['MV Name'])  
        MV_var.append('MV ' + str(data_MV[i]['MV Index']))

    # Other descriptions
    var_name, var_index = 2*CV_name + 2*MV_name, 2*CV_var + 2*MV_var
    limits = ['HL']*len(CV_name) + ['LL']*len(CV_name) + ['HL']*len(MV_name) + ['LL']*len(MV_name)

    # Insert into table
    table.insert(loc = 0, column = 'Limit', value = limits)
    table.insert(loc = 0, column = 'Name', value = var_name)
    table.insert(loc = 0, column = 'Variable', value = var_index)
    table = pd.DataFrame(data = table).sort_values(by = ['Perturbation', 'Lagrange', 'Relative error', 'Absolute error', 'Variable', 'Limit'], \
                                                    ascending = [False, False, False, False, True, True])
    priority = np.arange(1, len(var_name) + 1, dtype = int)
    table.insert(loc = 0, column = 'Priority', value = priority)

    return table


# Compare non-local analysis (perturbation solution) vs. Lagrange multipliers
def pert_analysis(pertR, pertI, rec, comp):
    try:
        data_types = ['NUMBER OF DATA', 'PERT (Perturbation) vs. REC (Relative loss)', 'PERT (Perturbation) vs. REC (Absolute loss)', \
                        'PERT (Lagrange) vs. REC (Relative loss)',  'PERT (Lagrange) vs. REC (Absolute loss)', 'PERT (Perturbation) vs. PERT (Lagrange)']
        pertR_N1 = pertR.sort_values(by = ['Perturbation', 'Variable', 'Limit'], ascending = [False, True, True])
        pertR_N2 = pertR.sort_values(by = ['Lagrange', 'Variable', 'Limit'], ascending = [False, True, True])
        pertI_N1 = pertI.sort_values(by = ['Perturbation', 'Variable', 'Limit'], ascending = [False, True, True])
        pertI_N2 = pertI.sort_values(by = ['Lagrange', 'Variable', 'Limit'], ascending = [False, True, True])
        rec_1 = rec.sort_values(by = ['Relative loss', 'Variable', 'Limit'], ascending = [False, True, True])
        rec_2 = rec.sort_values(by = ['Absolute loss', 'Variable', 'Limit'], ascending = [False, True, True])
        pert_N1 = [np.array(pertR_N1['Variable'] + pertR_N1['Limit']), np.array(pertI_N1['Variable'] + pertI_N1['Limit'])]
        pert_N2 = [np.array(pertR_N2['Variable'] + pertR_N2['Limit']), np.array(pertI_N2['Variable'] + pertI_N2['Limit'])]
        rec_N1 = np.array(rec_1['Variable'] + rec_1['Limit'])
        rec_N2 = np.array(rec_2['Variable'] + rec_2['Limit'])

        # Real model perturbation analysis
        Rmodel = [len(pertR), distance(pert_N1[0], rec_N1), distance(pert_N1[0], rec_N2), \
                    distance(pert_N2[0], rec_N1), distance(pert_N2[0], rec_N2), distance(pert_N1[0], pert_N2[0])]

        # Ideal model perturbation analysis
        Imodel = [len(pertI), distance(pert_N1[1], rec_N1), distance(pert_N1[1], rec_N2), \
                    distance(pert_N2[1], rec_N1), distance(pert_N2[1], rec_N2), distance(pert_N1[1], pert_N2[1])]

        # Table data
        table_ = pd.DataFrame(data = {'Edit distance': data_types, 'Real Model': Rmodel, 'Ideal Model': Imodel})
    except Exception as e:
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)

    return table_


# Optimal control point analysis
def table(case, args_data):

    try:
        # Load files
        fileExcel = re.search('SteadyStateBFCCU2C1.xlsx', args_data.xlsx).group(0)
        fileURT = re.search('BFCCU2C1_data(.+?).urt', args_data.urt).group(0)
        solution, solution2, solution3, solution4, M, M2, M3, M4, data_CV, data_MV, data = \
            load_from_file("data_" + fileURT + "_" + fileExcel + "_S.txt", 1)

    except Exception as e:
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)

    try:

        # Recommendations
        if case == 1:

            # Print saving data to file
            if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Recommendation is True:
                print(a[10])

            # Calculate recommendations
            table, table1 = check_limit(solution, solution2, solution3, solution4, M, M2, M3, M4, data_CV, data_MV, args_data)
            table, table1 = check_limit_comb(table, table1, data_CV, data_MV, solution, solution4, M2, M4)

            # Calculate constraint violations
            Table1 = check_violation(solution2, solution4, data_CV, data_MV)

            # Print data to terminal (console)
            if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Recommendation is True:
                print(), print(a[6])
                print(tabulate(Table1, headers = 'keys', floatfmt = ".4E", showindex = False))
                print(), print(a[8])
                print(tabulate(table, headers = 'keys', floatfmt = ".4E", showindex = False))
                print(), print(a[7])
                print(tabulate(table1, headers = 'keys', floatfmt = ".4E", showindex = False))
                print()

            # Save data to file
            if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Recommendation is True:

                table_name = 'Recommendations_' + fileURT.split('.')[0] + '_hard_S'
                table_name1 = 'Recommendations_' + fileURT.split('.')[0] + '_1_norm_S'

                # Save to .xlsx file
                if args_data.Format.lower() == 'xlsx':
                    df = {}
                    df.update(zip(["Legend", ""],[list(map(itemgetter(0), [[y[0].strip(), y[1].strip()] for y in [x.split("-") for x in a[1].split("\n")[3:-2]]])), \
                                                    list(map(itemgetter(1), [[y[0].strip(), y[1].strip()] for y in [x.split("-") for x in a[1].split("\n")[3:-2]]]))]))
                    df = pd.DataFrame(data = df)
                    with pd.ExcelWriter(table_name1 + ".xlsx") as writer1:
                        df.to_excel(writer1, sheet_name = 'Legend', index = False)
                        table1.to_excel(writer1, sheet_name = 'Recomms', index = False)
                        Table1.to_excel(writer1, sheet_name = 'Valids', index = False)
                    with pd.ExcelWriter(table_name + ".xlsx") as writer:
                        df.to_excel(writer, sheet_name = 'Legend', index = False)
                        table.to_excel(writer, sheet_name = 'Recomms', index = False)
                        Table1.to_excel(writer, sheet_name = 'Valids', index = False)

                # Save to .txt file
                elif args_data.Format.lower() == 'txt':
                    with open(table_name1 + ".txt", 'w') as file1:
                        file1.write(a[1])
                        file1.write(a[2])
                        file1.write(tabulate(Table1, headers = 'keys', showindex = False, tablefmt = args_data.Format.lower()))
                        file1.write(a[3])
                        file1.write(tabulate(table1, headers = 'keys', showindex = False, tablefmt = args_data.Format.lower()))
                        file1.close()
                    with open(table_name + ".txt", 'w') as file:
                        file.write(a[1])
                        file.write(a[2])
                        file.write(tabulate(Table1, headers = 'keys', showindex = False, tablefmt = args_data.Format.lower()))
                        file.write(a[3])
                        file.write(tabulate(table, headers = 'keys', showindex = False, tablefmt = args_data.Format.lower()))
                        file.close()

            # Open file/s
            if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Recommendation is True:
                if args_data.Format.lower() == 'xlsx':
                    os.system('start "" /max ' + table_name + '.xlsx')
                    os.system('start "" /max ' + table_name1 + '.xlsx')
                elif args_data.Format.lower() == 'txt':
                    os.system('start "" /max ' + table_name + '.txt')
                    os.system('start "" /max ' + table_name1 + '.txt')

            # Show info
            if args_data.Recommendation is True:
                if 'BFCCU2C1_data' in sys.argv[-2]:
                    print('\n!!!! Recommendations calculation solved !!!!\n\n')
                else:
                    print('\n' + stylize('Recommendations calculation solved', fg('green') + attr('bold')) + '\n\n')

    except Exception as e:
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)

    try:

        # Perturbations
        if case == 2:

            solvers.options['abstol'] = 1e-12
            solvers.options['reltol'] = 1e-12

            # Print saving data to file
            if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Perturbation is True:
                print(a[10])

            # Non-local analysis (perturbation solution) - change of the right side = 1e-4
            if args_data.Analysis == 'Real' and args_data.Perturbation is True:
                perturbations = perturbations_analysis(M['P'], M['q'], M['G'], M['A'], M['b'], solution['real'], 0, 0, data_CV, data_MV, args_data, solution4, M4)
                perturbations1 = perturbations_analysis(M4['P'], M4['q'], M4['G'], M4['A'], M4['b'], solution4['real'], 0, 1, data_CV, data_MV, args_data, solution4, M4)
                model_type = 'Real_'

            # Non-local analysis (perturbation solution) - change of the right side = abs(h_i,ideal - h_i,real)
            elif args_data.Analysis == 'Ideal' and args_data.Perturbation is True:
                perturbations = perturbations_analysis(M['P'], M['q'], M['G'], M['A'], M['b'], solution['real'], 1, 0, data_CV, data_MV, args_data, solution4, M4)
                perturbations1 = perturbations_analysis(M4['P'], M4['q'], M4['G'], M4['A'], M4['b'], solution4['real'], 1, 1, data_CV, data_MV, args_data, solution4, M4)
                model_type = 'Ideal_'

            # Non-local analysis (perturbation solution) - change of the right side = 1e-4 & change of the right side = abs(h_i,ideal - h_i,real)
            elif args_data.Analysis == 'Constraints' and args_data.Perturbation is True:
                perturbations_ = perturbations_analysis(M['P'], M['q'], M['G'], M['A'], M['b'], solution['real'], 0, 0, data_CV, data_MV, args_data, solution4, M4)
                _perturbations = perturbations_analysis(M['P'], M['q'], M['G'], M['A'], M['b'], solution['real'], 1, 0, data_CV, data_MV, args_data, solution4, M4)
                perturbations1_ = perturbations_analysis(M4['P'], M4['q'], M4['G'], M4['A'], M4['b'], solution4['real'], 0, 1, data_CV, data_MV, args_data, solution4, M4)
                _perturbations1 = perturbations_analysis(M4['P'], M4['q'], M4['G'], M4['A'], M4['b'], solution4['real'], 1, 1, data_CV, data_MV, args_data, solution4, M4)
                _table, _table1 = check_limit(solution, solution2, solution3, solution4, M, M2, M3, M4, data_CV, data_MV, args_data)
                model_type = 'Analysis_'

            # Non-local analysis (perturbation solution) - table data
            if (args_data.Analysis == 'Real' or args_data.Analysis == 'Ideal') and args_data.Perturbation is True:
                model = 1 if args_data.Analysis == 'Ideal' else 0
                table = table_perturbations(perturbations, data_CV, data_MV, model)
                table1 = table_perturbations(perturbations1, data_CV, data_MV, model)
            
            # Compare non-local analysis (perturbation solution) vs. Lagrange multipliers
            if args_data.Analysis == 'Constraints' and args_data.Perturbation is True:
                table__ = table_perturbations(perturbations_, data_CV, data_MV, 0)
                _table_ = table_perturbations(_perturbations, data_CV, data_MV, 1)
                table_ = pert_analysis(table__, _table_, _table, args_data.Rcomp)
                table_1_ = table_perturbations(perturbations1_, data_CV, data_MV, 0)
                _table_1 = table_perturbations(_perturbations1, data_CV, data_MV, 1)
                table1_ = pert_analysis(table_1_, _table_1, _table1, args_data.Rcomp)

            # Print data to terminal (console)
            if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Perturbation is True:
                if args_data.Analysis == 'Constraints':
                    print(), print(a[12])
                    print(tabulate(table_, headers = 'keys', showindex = False))
                    print(), print(a[11])
                    print(tabulate(table1_, headers = 'keys', showindex = False))
                    print()
                elif args_data.Analysis == 'Real':
                    print(), print(a[13])
                    print(tabulate(table, headers = 'keys', showindex = False))
                    print(), print(a[14])
                    print(tabulate(table1, headers = 'keys', showindex = False))
                    print()
                elif args_data.Analysis == 'Ideal':
                    print(), print(a[16])
                    print(tabulate(table, headers = 'keys', showindex = False))
                    print(), print(a[17])
                    print(tabulate(table1, headers = 'keys', showindex = False))
                    print()

            # Save data to file
            if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Perturbation is True:

                # Save to .xlsx file
                if args_data.Format.lower() == 'xlsx':
                    if args_data.Analysis == 'Real' or args_data.Analysis == 'Ideal':
                        a_n = 19
                    elif args_data.Analysis == 'Constraints':
                        a_n = 20
                    df = {}
                    df.update(zip(["Legend", ""],[list(map(itemgetter(0), [[y[0].strip(), y[1].strip()] for y in [x.split("-") for x in a[a_n].split("\n")[3:-2]]])), \
                                                    list(map(itemgetter(1), [[y[0].strip(), y[1].strip()] for y in [x.split("-") for x in a[a_n].split("\n")[3:-2]]]))]))
                    df = pd.DataFrame(data = df)
                    if (args_data.Analysis == 'Real' or args_data.Analysis == 'Ideal') and not type(perturbations) == type('string'):
                        table_name = 'Perturbations_' + model_type + fileURT.split('.')[0] + '_hard_S.xlsx'
                        with pd.ExcelWriter(table_name) as writer:
                            df.to_excel(writer, sheet_name = 'Legend', index = False)
                            table.to_excel(writer, sheet_name = 'Data', index = False)
                    if (args_data.Analysis == 'Real' or args_data.Analysis == 'Ideal') and not type(perturbations1) == type('string'):
                        table_name1 = 'Perturbations_' + model_type + fileURT.split('.')[0] + '_1_norm_S.xlsx'
                        with pd.ExcelWriter(table_name1) as writer1:
                            df.to_excel(writer1, sheet_name = 'Legend', index = False)
                            table1.to_excel(writer1, sheet_name = 'Data', index = False)
                    if args_data.Analysis == 'Constraints' and (not type(perturbations1_) == type('string') and not type(_perturbations1) == type('string')):
                        table_name1_ = 'Perturbations_' + model_type + fileURT.split('.')[0] + '_1_norm_S.xlsx'
                        with pd.ExcelWriter(table_name1_) as writer1_:
                            df.to_excel(writer1_, sheet_name = 'Legend', index = False)
                            table1_.to_excel(writer1_, sheet_name = 'Data', index = False)
                    if args_data.Analysis == 'Constraints' and (not type(perturbations_) == type('string') and not type(_perturbations) == type('string')):
                        table_name_ = 'Perturbations_' + model_type + fileURT.split('.')[0] + '_hard_S.xlsx'
                        with pd.ExcelWriter(table_name_) as writer_:
                            df.to_excel(writer_, sheet_name = 'Legend', index = False)
                            table_.to_excel(writer_, sheet_name = 'Data', index = False)

                # Save to .txt file
                elif args_data.Format.lower() == 'txt':
                    if (args_data.Analysis == 'Real' or args_data.Analysis == 'Ideal') and not type(perturbations) == type('string'):
                        table_name = 'Perturbations_' + model_type + fileURT.split('.')[0] + '_hard_S.txt'
                        with open(table_name, 'w') as file:
                            file.write(a[19] + '\n\n')
                            file.write(tabulate(table, headers = 'keys', showindex = False, tablefmt = args_data.Format.lower()))
                            file.close()
                    if (args_data.Analysis == 'Real' or args_data.Analysis == 'Ideal') and not type(perturbations1) == type('string'):
                        table_name1 = 'Perturbations_' + model_type + fileURT.split('.')[0] + '_1_norm_S.txt'
                        with open(table_name1, 'w') as file1:
                            file1.write(a[19] + '\n\n')
                            file1.write(tabulate(table1, headers = 'keys', showindex = False, tablefmt = args_data.Format.lower()))
                            file1.close()
                    if args_data.Analysis == 'Constraints' and (not type(perturbations1_) == type('string') and not type(_perturbations1) == type('string')):
                        table_name1_ = 'Perturbations_' + model_type + fileURT.split('.')[0] + '_1_norm_S.txt'
                        with open(table_name1_, 'w') as file1_:
                            file1_.write(a[19] + '\n\n')
                            file1_.write(tabulate(table1_, headers = 'keys', showindex = False, tablefmt = args_data.Format.lower()))
                            file1_.close()
                    if args_data.Analysis == 'Constraints' and (not type(perturbations_) == type('string') and not type(_perturbations) == type('string')):
                        table_name_ = 'Perturbations_' + model_type + fileURT.split('.')[0] + '_hard_S.txt'
                        with open(table_name_, 'w') as file_:
                            file_.write(a[19] + '\n\n')
                            file_.write(tabulate(table_, headers = 'keys', showindex = False, tablefmt = args_data.Format.lower()))
                            file_.close()

            # Open file/s
            if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Perturbation is True:
                if args_data.Analysis == 'Real' or args_data.Analysis == 'Ideal':
                    os.system('start "" /max ' + table_name)
                    os.system('start "" /max ' + table_name1)
                elif args_data.Analysis == 'Constraints':
                    os.system('start "" /max ' + table_name_)
                    os.system('start "" /max ' + table_name1_)

            # Show info
            if args_data.Perturbation is True:
                if args_data.Analysis == 'Ideal':
                    if 'BFCCU2C1_data' in sys.argv[-2]:
                        print('\n!!!! Ideal change perturbations calculation solved !!!!\n\n')
                    else:
                        print('\n' + stylize('Ideal change perturbations calculation solved', fg('green') + attr('bold')) + '\n\n')
                elif args_data.Analysis == 'Real':
                    if 'BFCCU2C1_data' in sys.argv[-2]:
                        print('\n!!!! Real change perturbations calculation solved !!!!\n\n')
                    else:
                        print('\n' + stylize('Real change perturbations calculation solved', fg('green') + attr('bold')) + '\n\n')
                elif args_data.Analysis == 'Constraints':
                    if 'BFCCU2C1_data' in sys.argv[-2]:
                        print('\n!!!! Constraints perturbations calculation solved !!!!\n\n')
                    else:
                        print('\n' + stylize('Constraints perturbations calculation solved', fg('green') + attr('bold')) + '\n\n')

    except Exception as e:
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)


# Optimal control point analysis
def table2(case, args_data):

    try:
        # Load files
        fileExcel = re.search('Zošit(.+?).xlsx', args_data.xlsx).group(0)
        fileURT = re.search('BFCCU2C1_data(.+?).urt', args_data.urt).group(0)
        config_values = pd.read_excel(fileExcel, 'Sheet1')
        solution, solution2, solution3, solution4, M, M2, M3, M4, data_CV, data_MV, data = \
            load_from_file("data_" + fileURT + "_" + fileExcel + "_M.txt", 1)

    except Exception as e:
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)

    try:

        # Recommendations
        if case == 1:

            # Print saving data to file
            if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Recommendation is True:
                print(a[10])

            # Iterate over time data
            table, table1, Table1 = {}, {}, {}
            for dtX in range(1, len(data.items()) + 1):

                # Set file open method
                h_var = "w" if dtX == 1 else "a"

                # Print current time period
                if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Recommendation is True:
                    print("\n" + a[5] + "\n")
                    print(config_values.columns[dtX].strftime("%Y-%m-%d %H:%M:%S"))
                    print("\n" + a[5] + "\n")

                # Calculate recommendations
                table[dtX], table1[dtX] = check_limit(solution[dtX], solution2[dtX], solution3[dtX], solution4[dtX], \
                                                M[dtX], M2[dtX], M3[dtX], M4[dtX], data_CV[dtX], data_MV[dtX], args_data)
                table[dtX], table1[dtX] = check_limit_comb(table[dtX], table1[dtX], data_CV[dtX], data_MV[dtX], solution[dtX], solution4[dtX], M[dtX], M4[dtX])

                # Calculate constraint violations
                Table1[dtX] = check_violation(solution2[dtX], solution4[dtX], data_CV[dtX], data_MV[dtX])

                # Print data to terminal (console)
                if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Recommendation is True:
                    print(), print(a[6])
                    print(tabulate(Table1[dtX][Table1[dtX]['Violation'] != 0], headers = 'keys', floatfmt = ".4E", showindex = False))
                    print(), print(a[8])
                    print(tabulate(table[dtX][table[dtX]['Relative loss'] != 0], headers = 'keys', floatfmt = ".4E", showindex = False))
                    print(), print(a[7])
                    print(tabulate(table1[dtX][table1[dtX]['Relative loss'] != 0], headers = 'keys', floatfmt = ".4E", showindex = False))
                    print()

                # Save data to file
                if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Recommendation is True:

                    table_name = 'Recommendations_' + fileURT.split('.')[0] + '_hard_M'
                    table_name1 = 'Recommendations_' + fileURT.split('.')[0] + '_1_norm_M'

                    # Save to .xlsx file
                    if args_data.Format.lower() == 'xlsx':
                        pass
                        if dtX == 1:
                            df = {}
                            df.update(zip(["Legend", ""],[list(map(itemgetter(0), [[y[0].strip(), y[1].strip()] for y in [x.split("-") for x in a[1].split("\n")[3:-2]]])), \
                                                         list(map(itemgetter(1), [[y[0].strip(), y[1].strip()] for y in [x.split("-") for x in a[1].split("\n")[3:-2]]]))]))
                            df = pd.DataFrame(data = df)
                        with pd.ExcelWriter(table_name + ".xlsx", engine = "openpyxl", mode = h_var) as writer:
                            if dtX == 1: df.to_excel(writer, sheet_name = 'Legend', index = False)
                            table[dtX][table[dtX]['Relative loss'] > 0].to_excel(writer, sheet_name = 'Recomms_' + str(config_values.columns[dtX].strftime("%Y.%m.%d %H.%M.%S")), index = False)
                            Table1[dtX][Table1[dtX]['Violation'] > 0].to_excel(writer, sheet_name = 'Valids_' + str(config_values.columns[dtX].strftime("%Y.%m.%d %H.%M.%S")), index = False)
                        with pd.ExcelWriter(table_name1 + ".xlsx", engine = "openpyxl", mode = h_var) as writer1:
                            if dtX == 1: df.to_excel(writer1, sheet_name = 'Legend', index = False)
                            table1[dtX][table1[dtX]['Relative loss'] > 0].to_excel(writer1, sheet_name = 'Recomms_' + str(config_values.columns[dtX].strftime("%Y.%m.%d %H.%M.%S")), index = False)
                            Table1[dtX][Table1[dtX]['Violation'] > 0].to_excel(writer, sheet_name = 'Valids_' + str(config_values.columns[dtX].strftime("%Y.%m.%d %H.%M.%S")), index = False)

                    # Save to .txt file
                    elif args_data.Format.lower() == 'txt':
                        with open(table_name + ".txt", h_var) as file:
                            if dtX == 1: file.write(a[1])
                            if dtX != 1: file.write('\n')
                            file.write('\n\n' + a[4] + str(config_values.columns[dtX].strftime("%Y-%m-%d %H:%M:%S")) + a[4])
                            file.write(a[2])
                            file.write(tabulate(Table1[dtX][Table1[dtX]['Violation'] > 0], headers = 'keys', showindex = False, tablefmt = args_data.Format.lower()))
                            file.write(a[3])
                            file.write(tabulate(table[dtX][table[dtX]['Relative loss'] > 0], headers = 'keys', showindex = False, tablefmt = args_data.Format.lower()))
                            file.close()
                        with open(table_name1 + ".txt", h_var) as file1:
                            if dtX == 1: file1.write(a[1])
                            if dtX != 1: file1.write('\n')
                            file1.write('\n\n' + a[4] + str(config_values.columns[dtX].strftime("%Y-%m-%d %H:%M:%S")) + a[4])
                            file1.write(a[2])
                            file1.write(tabulate(Table1[dtX][Table1[dtX]['Violation'] > 0], headers = 'keys', showindex = False, tablefmt = args_data.Format.lower()))
                            file1.write(a[3])
                            file1.write(tabulate(table1[dtX][table1[dtX]['Relative loss'] > 0], headers = 'keys', showindex = False, tablefmt = args_data.Format.lower()))
                            file1.close()

            # Open file/s
            if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Recommendation is True:
                if args_data.Format.lower() == 'xlsx':
                    os.system('start "" /max ' + table_name + '.xlsx')
                    os.system('start "" /max ' + table_name1 + '.xlsx')
                elif args_data.Format.lower() == 'txt':
                    os.system('start "" /max ' + table_name + '.txt')
                    os.system('start "" /max ' + table_name1 + '.txt')

            # Show info
            if args_data.Recommendation is True:
                if 'BFCCU2C1_data' in sys.argv[-2]:
                    print('\n!!!! Recommendations calculation solved !!!!\n\n')
                else:
                    print('\n' + stylize('Recommendations calculation solved', fg('green') + attr('bold')) + '\n\n')
            if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Recommendation is True:
                if 'BFCCU2C1_data' in sys.argv[-2]:
                    print('\n!!!! Data was saved to external file !!!!\n\n')
                else:
                    print('\n' + stylize('Data was saved to external file', fg('blue') + attr('bold')) + '\n\n')

            # Plot overall loss data
            datasA, datasR, names = {}, {}, sorted(list(table[1]['Variable']+' '+table[1]['Limit'])[1:])
            for i in names:
                datasA[i] = [0.0]*len(data.items())
                datasR[i] = [0.0]*len(data.items())
            for i in range(1, len(data.items()) + 1):
                for _, j in table[i].iterrows():
                    if j['Variable'] + ' ' + j['Limit'] == "   ": continue
                    if j['Absolute loss'] > 0.0:
                        datasA[j['Variable'] + ' ' + j['Limit']][i-1] = j['Absolute loss']
                    if j['Relative loss'] > 0.0:
                        datasR[j['Variable'] + ' ' + j['Limit']][i-1] = j['Relative loss']

            vrb, limit, cnt, initial, final = [], [], [], [], []
            for key, value in datasA.items():
                idx = [i for i, e in enumerate(value) if e != 0.0]
                try:
                    t_init, t_final, counter = idx[0], idx[0], 1
                except:
                    t_init, t_final, counter = 0, 0, 1
                for pos in range(len(idx)-1):
                    if idx[pos+1] - idx[pos] <= 3:
                        t_final = idx[pos+1]
                    else:
                        if t_final - t_init >= args_data.Interval:
                            vrb.append(key.split()[0] + " " + key.split()[1])
                            limit.append(key.split()[2])
                            cnt.append(counter)
                            initial.append(t_init)
                            final.append(t_final)
                            counter += 1
                        t_init, t_final = idx[pos+1], idx[pos+1]
                    if pos == len(idx)-2:
                        if t_final - t_init >= args_data.Interval:
                            vrb.append(key.split()[0] + " " + key.split()[1])
                            limit.append(key.split()[2])
                            cnt.append(counter)
                            initial.append(t_init)
                            final.append(t_final)
                            counter += 1
            time = [config_values.columns[i] for i in range(1, len(data.items()) + 1)]
            problems = {"Variable": vrb, "Limit": limit, "Interval": cnt, "Initial Time": [time[i] for i in initial], "Final Time": [time[i] for i in final]}
            print(tabulate(pd.DataFrame(data = problems), headers = 'keys', showindex = False))

            plt.rcParams['text.usetex'] = True
            plt.figure(1)
            for key, value in datasA.items():
                if not all(v == 0.0 for v in value):
                    plt.step([config_values.columns[i] for i in range(1, len(data.items()) + 1)], value, label = key)
            plt.subplots_adjust(bottom = 0.25)
            plt.grid(linewidth = 0.5)
            plt.yscale('log')
            plt.ylabel("Absolute loss")
            plt.legend()
            plt.figure(2)
            for key, value in datasR.items():
                if not all(v == 0.0 for v in value):
                    plt.step([config_values.columns[i] for i in range(1, len(data.items()) + 1)], value, label = key)
            plt.subplots_adjust(bottom = 0.25)
            plt.grid(linewidth = 0.5)
            plt.yscale('log')
            plt.ylabel("Relative loss (in \%)")
            plt.legend()
            plt.show()

    except Exception as e:
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)

    try:

        # Perturbations
        if case == 2:

            solvers.options['abstol'] = 1e-12
            solvers.options['reltol'] = 1e-12

            # Print saving data to file
            if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Perturbation is True:
                print(a[10]), print()

            deg_free = [0.0]*len(data.items())
            for dtX in range(1, len(data.items()) + 1):

                # Set file open method
                h_var = "w" if dtX == 1 else "a"

                # Print current time period
                if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Perturbation is True:
                    print("\n" + a[5] + "\n")
                    print(config_values.columns[dtX].strftime("%Y-%m-%d %H:%M:%S"))
                    print("\n" + a[5] + "\n")

                # Non-local analysis (perturbation solution) - change of the right side = 1e-4
                if args_data.Analysis == 'Real' and args_data.Perturbation is True:
                    perturbations = perturbations_analysis(M[dtX]['P'], M[dtX]['q'], M[dtX]['G'], M[dtX]['A'], M[dtX]['b'], solution[dtX]['real'], 0, 0, data_CV[dtX], data_MV[dtX], args_data, solution4[dtX], M4[dtX])
                    perturbations1 = perturbations_analysis(M4[dtX]['P'], M4[dtX]['q'], M4[dtX]['G'], M4[dtX]['A'], M4[dtX]['b'], solution4[dtX]['real'], 0, 1, data_CV[dtX], data_MV[dtX], args_data, solution4[dtX], M4[dtX])
                    model_type = 'Real_'

                # Non-local analysis (perturbation solution) - change of the right side = abs(h_i,ideal - h_i,real)
                elif args_data.Analysis == 'Ideal' and args_data.Perturbation is True:
                    perturbations = perturbations_analysis(M[dtX]['P'], M[dtX]['q'], M[dtX]['G'], M[dtX]['A'], M[dtX]['b'], solution[dtX]['real'], 1, 0, data_CV[dtX], data_MV[dtX], args_data, solution4[dtX], M4[dtX])
                    perturbations1 = perturbations_analysis(M4[dtX]['P'], M4[dtX]['q'], M4[dtX]['G'], M4[dtX]['A'], M4[dtX]['b'], solution4[dtX]['real'], 1, 1, data_CV[dtX], data_MV[dtX], args_data, solution4[dtX], M4[dtX])
                    model_type = 'Ideal_'

                # Non-local analysis (perturbation solution) - change of the right side = 1e-4 & change of the right side = abs(h_i,ideal - h_i,real)
                elif args_data.Analysis == 'Constraints' and args_data.Perturbation is True:
                    perturbations_ = perturbations_analysis(M[dtX]['P'], M[dtX]['q'], M[dtX]['G'], M[dtX]['A'], M[dtX]['b'], solution[dtX]['real'], 0, 0, data_CV[dtX], data_MV[dtX], args_data, solution4[dtX], M4[dtX])
                    _perturbations = perturbations_analysis(M[dtX]['P'], M[dtX]['q'], M[dtX]['G'], M[dtX]['A'], M[dtX]['b'], solution[dtX]['real'], 1, 0, data_CV[dtX], data_MV[dtX], args_data, solution4[dtX], M4[dtX])
                    perturbations1_ = perturbations_analysis(M4[dtX]['P'], M4[dtX]['q'], M4[dtX]['G'], M4[dtX]['A'], M4[dtX]['b'], solution4[dtX]['real'], 0, 1, data_CV[dtX], data_MV[dtX], args_data, solution4[dtX], M4[dtX])
                    _perturbations1 = perturbations_analysis(M4[dtX]['P'], M4[dtX]['q'], M4[dtX]['G'], M4[dtX]['A'], M4[dtX]['b'], solution4[dtX]['real'], 1, 1, data_CV[dtX], data_MV[dtX], args_data, solution4[dtX], M4[dtX])
                    _table, _table1 = check_limit(solution[dtX], solution2[dtX], solution3[dtX], solution4[dtX], M[dtX], M2[dtX], M3[dtX], M4[dtX], data_CV[dtX], data_MV[dtX], args_data)
                    model_type = 'Analysis_'

                # Non-local analysis (perturbation solution) - table data
                if (args_data.Analysis == 'Real' or args_data.Analysis == 'Ideal') and args_data.Perturbation is True:
                    model = 1 if args_data.Analysis == 'Ideal' else 0
                    table = table_perturbations(perturbations, data_CV[dtX], data_MV[dtX], model)
                    table1 = table_perturbations(perturbations1, data_CV[dtX], data_MV[dtX], model)
                
                # Compare non-local analysis (perturbation solution) vs. Lagrange multipliers
                if args_data.Analysis == 'Constraints' and args_data.Perturbation is True:
                    table__ = table_perturbations(perturbations_, data_CV[dtX], data_MV[dtX], 0)
                    _table_ = table_perturbations(_perturbations, data_CV[dtX], data_MV[dtX], 1)
                    table_ = pert_analysis(table__, _table_, _table, args_data.Rcomp)
                    table_1_ = table_perturbations(perturbations1_, data_CV[dtX], data_MV[dtX], 0)
                    _table_1 = table_perturbations(_perturbations1, data_CV[dtX], data_MV[dtX], 1)
                    table1_ = pert_analysis(table_1_, _table_1, _table1, args_data.Rcomp)

                # Print data to terminal (console)
                if (args_data.Display == 'Terminal' or args_data.Display == 'Both') and args_data.Perturbation is True:
                    if args_data.Analysis == 'Constraints':
                        print(), print(a[12])
                        print(tabulate(table_, headers = 'keys', showindex = False))
                        print(), print(a[11])
                        print(tabulate(table1_, headers = 'keys', showindex = False))
                        print()
                    elif args_data.Analysis == 'Real':
                        print(), print(a[13])
                        print(tabulate(table, headers = 'keys', showindex = False))
                        print(), print(a[14])
                        print(tabulate(table1, headers = 'keys', showindex = False))
                        print()
                    elif args_data.Analysis == 'Ideal':
                        print(), print(a[16])
                        print(tabulate(table, headers = 'keys', showindex = False))
                        print(), print(a[17])
                        print(tabulate(table1, headers = 'keys', showindex = False))
                        print()

                # Save data to file
                if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Perturbation is True:

                    # Save to .xlsx file
                    if args_data.Format.lower() == 'xlsx':
                        if dtX == 1:
                            if args_data.Analysis == 'Real' or args_data.Analysis == 'Ideal':
                                a_n = 19
                            elif args_data.Analysis == 'Constraints':
                                a_n = 20
                            df = {}
                            df.update(zip(["Legend", ""],[list(map(itemgetter(0), [[y[0].strip(), y[1].strip()] for y in [x.split("-") for x in a[a_n].split("\n")[3:-2]]])), \
                                                         list(map(itemgetter(1), [[y[0].strip(), y[1].strip()] for y in [x.split("-") for x in a[a_n].split("\n")[3:-2]]]))]))
                            df = pd.DataFrame(data = df)
                        if (args_data.Analysis == 'Real' or args_data.Analysis == 'Ideal') and not type(perturbations) == type('string'):
                            table_name = 'Perturbations_' + model_type + fileURT.split('.')[0] + '_M.xlsx'
                            with pd.ExcelWriter(table_name + ".xlsx", engine = "openpyxl", mode = h_var) as writer:
                                if dtX == 1:
                                    df.to_excel(writer, sheet_name = 'Legend', index = False)
                                table.to_excel(writer, sheet_name = 'Data_' + str(config_values.columns[dtX].strftime("%Y.%m.%d %H.%M.%S")), index = False)
                        if (args_data.Analysis == 'Real' or args_data.Analysis == 'Ideal') and not type(perturbations1) == type('string'):
                            table_name1 = 'Perturbations_' + model_type + fileURT.split('.')[0] + '_1_norm_M.xlsx'
                            with pd.ExcelWriter(table_name1 + ".xlsx", engine = "openpyxl", mode = h_var) as writer1:
                                if dtX == 1:
                                    df.to_excel(writer1, sheet_name = 'Legend', index = False)
                                table1.to_excel(writer1, sheet_name = 'Data_' + str(config_values.columns[dtX].strftime("%Y.%m.%d %H.%M.%S")), index = False)
                        if args_data.Analysis == 'Constraints' and (not type(perturbations1_) == type('string') and not type(_perturbations1) == type('string')):
                            table_name1_ = 'Perturbations_' + model_type + fileURT.split('.')[0] + '_1_norm_M.xlsx'
                            with pd.ExcelWriter(table_name1_ + ".xlsx", engine = "openpyxl", mode = h_var) as writer1_:
                                if dtX == 1:
                                    df.to_excel(writer1_, sheet_name = 'Legend', index = False)
                                table1_.to_excel(writer1_, sheet_name = 'Data_' + str(config_values.columns[dtX].strftime("%Y.%m.%d %H.%M.%S")), index = False)
                        if args_data.Analysis == 'Constraints' and (not type(perturbations_) == type('string') and not type(_perturbations) == type('string')):
                            table_name_ = 'Perturbations_' + model_type + fileURT.split('.')[0] + '_M.xlsx'
                            with pd.ExcelWriter(table_name_ + ".xlsx", engine = "openpyxl", mode = h_var) as writer_:
                                if dtX == 1:
                                    df.to_excel(writer_, sheet_name = 'Legend', index = False)
                                table_.to_excel(writer_, sheet_name = 'Data_' + str(config_values.columns[dtX].strftime("%Y.%m.%d %H.%M.%S")), index = False)

                    # Save to .txt file
                    elif args_data.Format.lower() == 'txt':
                        if (args_data.Analysis == 'Real' or args_data.Analysis == 'Ideal') and not type(perturbations) == type('string'):
                            table_name = 'Perturbations_' + model_type + fileURT.split('.')[0] + '_M.txt'
                            with open(table_name, h_var) as file:
                                if dtX == 1: file.write(a[19])
                                if dtX != 1: file.write('\n')
                                file.write('\n\n' + a[4] + str(config_values.columns[dtX].strftime("%Y-%m-%d %H:%M:%S")) + a[4] + '\n')
                                file.write(tabulate(table, headers = 'keys', showindex = False, tablefmt = args_data.Format.lower()))
                                file.close()
                        if (args_data.Analysis == 'Real' or args_data.Analysis == 'Ideal') and not type(perturbations1) == type('string'):
                            table_name1 = 'Perturbations_' + model_type + fileURT.split('.')[0] + '_1_norm_M.txt'
                            with open(table_name1, h_var) as file1:
                                if dtX == 1: file1.write(a[19])
                                if dtX != 1: file1.write('\n')
                                file1.write('\n\n' + a[4] + str(config_values.columns[dtX].strftime("%Y-%m-%d %H:%M:%S")) + a[4] + '\n')
                                file1.write(tabulate(table1, headers = 'keys', showindex = False, tablefmt = args_data.Format.lower()))
                                file1.close()
                        if args_data.Analysis == 'Constraints' and (not type(perturbations1_) == type('string') and not type(_perturbations1) == type('string')):
                            table_name1_ = 'Perturbations_' + model_type + fileURT.split('.')[0] + '_1_norm_M.txt'
                            with open(table_name1_, h_var) as file1_:
                                if dtX == 1: file1_.write(a[20])
                                if dtX != 1: file1_.write('\n')
                                file1_.write('\n\n' + a[4] + str(config_values.columns[dtX].strftime("%Y-%m-%d %H:%M:%S")) + a[4] + '\n')
                                file1_.write(tabulate(table1_, headers = 'keys', showindex = False, tablefmt = args_data.Format.lower()))
                                file1_.close()
                        if args_data.Analysis == 'Constraints' and (not type(perturbations_) == type('string') and not type(_perturbations) == type('string')):
                            table_name_ = 'Perturbations_' + model_type + fileURT.split('.')[0] + '_M.txt'
                            with open(table_name_, h_var) as file_:
                                if dtX == 1: file_.write(a[20])
                                if dtX != 1: file_.write('\n')
                                file_.write('\n\n' + a[4] + str(config_values.columns[dtX].strftime("%Y-%m-%d %H:%M:%S")) + a[4] + '\n')
                                file_.write(tabulate(table_, headers = 'keys', showindex = False, tablefmt = args_data.Format.lower()))
                                file_.close()
                DF_help = table[table['Perturbation'] > 0.0]
                deg_free[dtX-1] = len(data_MV[dtX]) - len(DF_help[DF_help['Variable'].str.contains("MV")]) - len(DF_help[DF_help['Variable'].str.contains("CV")])
                break

            # Open file/s
            if (args_data.Display == 'File' or args_data.Display == 'Both') and args_data.Perturbation is True:
                if args_data.Analysis == 'Real' or args_data.Analysis == 'Ideal':
                    os.system('start "" /max ' + table_name)
                    os.system('start "" /max ' + table_name1)
                elif args_data.Analysis == 'Constraints':
                    os.system('start "" /max ' + table_name1_)
                    os.system('start "" /max ' + table_name_)

            # Show info
            if args_data.Perturbation is True:
                if args_data.Analysis == 'Ideal':
                    if 'BFCCU2C1_data' in sys.argv[-2]:
                        print('\n!!!! Ideal change perturbations calculation solved !!!!\n\n')
                    else:
                        print('\n' + stylize('Ideal change perturbations calculation solved', fg('green') + attr('bold')) + '\n\n')
                if args_data.Analysis == 'Real':
                    if 'BFCCU2C1_data' in sys.argv[-2]:
                        print('\n!!!! Real change perturbations calculation solved !!!!\n\n')
                    else:
                        print('\n' + stylize('Real change perturbations calculation solved', fg('green') + attr('bold')) + '\n\n')
                if args_data.Analysis == 'Constraints':
                    if 'BFCCU2C1_data' in sys.argv[-2]:
                        print('\n!!!! Constraints perturbations calculation solved !!!!\n\n')
                    else:
                        print('\n' + stylize('Constraints perturbations calculation solved', fg('green') + attr('bold')) + '\n\n')

            # Plot degrees of freedom
            deg_free = np.array(deg_free)
            deg_free2 = np.append(deg_free[1:],deg_free[-1])
            deg_diff = deg_free2 - deg_free
            plt.figure(3)
            plt.step([config_values.columns[i] for i in range(1, len(data.items()) + 1)], deg_free, '--o', color = 'black')
            plt.grid(linewidth = 0.5)
            plt.ylabel("Degrees of freedom $N_{DF}$", fontsize=12)
            plt.figure(4)
            from datetime import datetime
            plt.plot([[config_values.columns[i] for i in range(1, len(data.items()) + 1)][index] for index in [i for i, x in enumerate(deg_diff) if x > 0]], [x for x in deg_diff if x > 0], 'o', color = 'tab:green')
            plt.plot([[config_values.columns[i] for i in range(1, len(data.items()) + 1)][index] for index in [i for i, x in enumerate(deg_diff) if x < 0]], [x for x in deg_diff if x < 0], 'o', color = 'tab:red')
            plt.plot([[config_values.columns[i] for i in range(1, len(data.items()) + 1)][index] for index in [i for i, x in enumerate(deg_diff) if x == 0]], [x for x in deg_diff if x == 0], 'o', color = 'tab:blue')
            plt.grid(linewidth = 0.5)
            plt.ylabel("Change of degree of freedom d$N_{DF}$", fontsize=12)

            plt.show()

    except Exception as e:
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)


@Gooey(program_name="ORG - GUI/CLI", clear_before_run=True, terminal_font_family='Courier New')

def get_args():

    parser = GooeyParser(description='Tool for analyzing the reachability of the optimal operating point in APC programmed in the Python environment.')

    files = parser.add_argument_group('Load files',
                                      gooey_options={
                                        'show_border': True
    })

    files.add_argument("urt", widget="FileChooser",
                        help="Choose URT file",
                            gooey_options={
                                'wildcard':
                                    "URT files (*.urt)|*.urt|"
                                    "All files (*.*)|*.*"
                            }
    )

    files.add_argument("xlsx", widget="FileChooser",
                        help="Choose Excel file",
                            gooey_options={
                                'wildcard':
                                    "XLSX files (*.xlsx)|*.xlsx|"
                                    "All files (*.*)|*.*"
                            }
    )

    manipulation = parser.add_argument_group('Process data',
                                      gooey_options={
                                        'show_border': True
    })

    Todo = manipulation.add_mutually_exclusive_group()

    Todo.add_argument('-s', '--solve',
                      help="Find optimal control point",
                      action='store_true',
                      dest='Solve')

    Todo.add_argument('-r', '--recommendation',
                      help="Data validation and recommendations",
                      action='store_true',
                      dest='Recommendation')

    Todo.add_argument('-p', '--perturbation',
                      help="Model perturbations",
                      action='store_true',
                      dest='Perturbation')
    
    Poptions = parser.add_argument_group('Program settings',
                                      gooey_options={
                                        'show_border': True
    })

    Poptions.add_argument("-t", "--time",
                        choices=["Single", "Multiple"],
                        default="Single",
                        help="Number of data sets (default: Single)",
                        dest='Time')

    Poptions.add_argument("-a", "--analysis",
                        choices=["Real", "Ideal", "Constraints"],
                        default="Real",
                        help="Perturbation analysis (default: Real)",
                        dest='Analysis')

    Poptions.add_argument("-d", "--display",
                        choices=["Terminal", "File", "Both"],
                        default="Terminal",
                        help="Display data (default: Terminal)",
                        dest='Display')

    Poptions.add_argument("-f", "--format",
                        choices=["Txt", "Xlsx"],
                        default="Txt",
                        help="File format (default: Text)",
                        dest='Format')

    Poptions.add_argument("-g", "--graph",
                        action="store_true",
                        default=True,
                        help="Display graphs (default: True)",
                        dest='Graph')

    Poptions.add_argument("-ac", "--acomp",
                        type=float,
                        default=100.0,
                        help="Perturbation absolute comparison (default: 100.0)",
                        dest='Acomp')

    Poptions.add_argument("-rc", "--rcomp",
                        type=float,
                        default=0.01,
                        help="Perturbation relative comparison (default: 0.01)",
                        dest='Rcomp')
    
    Poptions.add_argument("-in", "--interval",
                        type=int,
                        default=10,
                        help="Constraint loss interval (default: 10)",
                        dest='Interval')

    Poptions.add_argument("-ch", "--change",
                        type=float,
                        default=1e-4,
                        help="Perturbation change (default: 1e-4)",
                        dest='Change')

    Poptions.add_argument("-ra", "--range",
                        type=float,
                        default=5.0,
                        help="Constraints range (default: 5.0)",
                        dest='Range')

    Soptions = parser.add_argument_group('Solver settings',
                                      gooey_options={
                                        'show_border': True
    })

    Soptions.add_argument("-at", "--atol",
                        type=float,
                        default=1e-7,
                        help="Absolute tolerance (default: 1e-7)",
                        dest='Atol')
    
    Soptions.add_argument("-rt", "--rtol",
                        type=float,
                        default=1e-6,
                        help="Relative tolerance (default: 1e-6)",
                        dest='Rtol')
    
    Soptions.add_argument("-mi", "--maxiter",
                        type=int,
                        default=100,
                        help="Maximal iterations (default: 100)",
                        dest='Maxiter')

    Soptions.add_argument("-pr", "--print",
                        action="store_true",
                        default=False,
                        help="Print solving iterations (default: False)",
                        dest='Print')

    return parser.parse_args()


def main():
    global a, t1, t2, t3
    with open("text.txt") as my_file:
        a = (my_file.read().split('#'))

    args = get_args()

    solvers.options['abstol'] = args.Atol
    solvers.options['reltol'] = args.Rtol
    solvers.options['maxiters'] = args.Maxiter
    solvers.options['show_progress'] = args.Print
    solvers.options['mosek'] = {iparam.log: int(args.Print), iparam.max_num_warnings: int(args.Print)}

    try:
        data_urt = re.search('BFCCU2C1_data(.+?).urt', args.urt).group(0)
    except:
        if 'BFCCU2C1_data' in sys.argv[-2]:
            print('\n!!!! Incorrectly entered/selected URT file !!!!\n\n')
        else:
            print('\n' + stylize('Incorrectly entered/selected URT file', fg('red') + attr('bold')) + '\n\n')
    try:
        if args.Time == 'Single' and (args.Solve is True or args.Recommendation is True or args.Perturbation is True):
            data_xlsx = re.search('SteadyStateBFCCU2C1.xlsx', args.xlsx).group(0)
        elif args.Time == 'Multiple' and (args.Solve is True or args.Recommendation is True or args.Perturbation is True):
            data_xlsx = re.search('Zošit(.+?).xlsx', args.xlsx).group(0)
        else:
            data_xlsx = args.xlsx.split("\\")[-1]
    except:
        if 'BFCCU2C1_data' in sys.argv[-2]:
            print('\n!!!! Incorrectly entered/selected XLSX file !!!!\n\n')
            exit()
        else:
            print('\n' + stylize('Incorrectly entered/selected XLSX file', fg('red') + attr('bold')) + '\n\n')
            exit()

    if args.Solve is True:
        option_TODO = "Solve"
    elif args.Recommendation is True:
        option_TODO = "Recommendation"
    elif args.Perturbation is True:
        option_TODO = "Perturbation"
    else:
        option_TODO = "None"
    print("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
    print("Set parameters\n==============")

    data_info = {"": ["Excel file", "URT file", "Choosen option", "Time", "Analysis", "Display", "Format", "Acomp", "Rcomp","Interval", "Change", "Range", \
                      "Atol", "Rtol", "Maxiter", "Print"], " ": [data_xlsx, data_urt, option_TODO, args.Time, args.Analysis, args.Display, args.Format, \
                      str(args.Acomp), str(args.Rcomp), str(args.Interval), str(args.Change), str(args.Range), str(args.Atol), str(args.Rtol), str(args.Maxiter), str(args.Print)]}
    print(tabulate(data_info, headers = 'keys'))
    print("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")

    if args.Time == 'Single':
        if args.Solve is True:
            execute(args)
        elif args.Recommendation is True:
            if os.path.isfile("data_" + data_urt + "_" + data_xlsx + "_S.txt") is True:
                t1, t2 = load_from_file("data_" + data_urt + "_" + data_xlsx + "_S.txt", 0)
                (_, _, _, _, _, _, _, _, mtime_URT, _) = os.stat(data_urt)
                (_, _, _, _, _, _, _, _, mtime_XLSX, _) = os.stat(data_xlsx)
                if t1 < int(mtime_URT) or t2 < int(mtime_XLSX):
                    execute(args)
            else:
                execute(args)
            table(1, args)
        elif args.Perturbation is True:
            execute(args), table(2, args)
    elif args.Time == 'Multiple':
        if args.Solve is True:
            execute2(args)
        elif args.Recommendation is True:
            if os.path.isfile("data_" + data_urt + "_" + data_xlsx + "_M.txt") is True:
                t1, t2 = load_from_file("data_" + data_urt + "_" + data_xlsx + "_M.txt", 0)
                (_, _, _, _, _, _, _, _, mtime_URT, _) = os.stat(data_urt)
                (_, _, _, _, _, _, _, _, mtime_XLSX, _) = os.stat(data_xlsx)
                if t1 < int(mtime_URT) or t2 < int(mtime_XLSX):
                    execute2(args)
            else:
                execute2(args)
            table2(1, args)
        elif args.Perturbation is True:
            if os.path.isfile("data_" + data_urt + "_" + data_xlsx + "_M.txt") is True:
                t1, t2 = load_from_file("data_" + data_urt + "_" + data_xlsx + "_M.txt", 0)
                (_, _, _, _, _, _, _, _, mtime_URT, _) = os.stat(data_urt)
                (_, _, _, _, _, _, _, _, mtime_XLSX, _) = os.stat(data_xlsx)
                if t1 < int(mtime_URT) or t2 < int(mtime_XLSX):
                    execute2(args)
            else:
                execute2(args)
            table2(2, args)

if __name__ == "__main__":
    main()