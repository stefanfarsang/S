import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


currentfolder = '/Users/stefanfarsang/Desktop/Test'  # folder for saving data, change it accordingly


p_bar = 2000  # pressure in barsgit status
t_c = 875  # temperature in Celsius (°C)
t_k = t_c + 273.15  # temperature in Kelvin (K)
fH2O = 1600.41  # H2O fugacity in bar from https://fluid-eos.web.psi.ch/EOS/calculator_simple.html
fc_H2S = 1  # fugacity coefficient of H2S
fc_SO2 = 1  # fugacity coefficient of SO2


# A function to calculate fO2 corresponding to the FMQ buffer
def fO2_fmq(t_k, p_bar):
        A = -25096.3
        B = 8.735
        C = 0.11
        logfO2_fmq = A/t_k + B + C * (p_bar-1)/t_k
        return logfO2_fmq


# A function to calculate fO2 corresponding to the NNO buffer
def fO2_nno(t_k, p_bar):
        A = -24930
        B = 9.36
        C = 0.046
        logfO2_nno = A/t_k + B + C * (p_bar-1)/t_k
        return logfO2_nno


# A function to calculate the SO2/H2S ratio
def fluid_s_speciation(t_k, fO2, fH2O, fc_SO2, fc_H2S):
        logK = +4.13 + 19516/t_k # from Figure S4c
        K = 10 ** logK
        SO2_to_H2S_ratio = K * (fO2 ** 1.5) / fH2O * (fc_H2S / fc_SO2)
        SO2_to_H2S_ratio = SO2_to_H2S_ratio / (SO2_to_H2S_ratio + 1)
        return SO2_to_H2S_ratio


# A function to calculate the speciation curve
def fluid_s_speciation_curve(buffer):
        oxygenfugacity = []
        sulfide_to_total_S = []
        SO2_to_total_S = []
        logfO2_from = -4
        logfO2_to = +4
        if buffer == 'FMQ':
                for i in range(logfO2_from * 1000, logfO2_to * 1000, 1):
                        oxygenfugacity.append(i / 1000)
                        result = fO2_fmq(t_k, p_bar)
                        logfO2 = result + i / 1000
                        SO2_to_total_S_per_point = fluid_s_speciation(t_k, 10 ** logfO2, fH2O, fc_SO2, fc_H2S)
                        SO2_to_total_S.append(SO2_to_total_S_per_point)
                        sulfide_to_total_S_per_point = 1 - fluid_s_speciation(t_k, 10 ** logfO2, fH2O, fc_SO2, fc_H2S)
                        sulfide_to_total_S.append(sulfide_to_total_S_per_point)
        if buffer == 'NNO':
                for i in range(logfO2_from*1000, logfO2_to*1000, 1):
                        oxygenfugacity.append(i/1000)
                        result = fO2_nno(t_k, p_bar)
                        logfO2 = result + i / 1000
                        SO2_to_total_S_per_point = fluid_s_speciation(t_k, 10 ** logfO2, fH2O, fc_SO2, fc_H2S)
                        SO2_to_total_S.append(SO2_to_total_S_per_point)
                        sulfide_to_total_S_per_point = 1 - fluid_s_speciation(t_k, 10 ** logfO2, fH2O, fc_SO2, fc_H2S)
                        sulfide_to_total_S.append(sulfide_to_total_S_per_point)
        return oxygenfugacity, sulfide_to_total_S, SO2_to_total_S


# A function to save figure as png, tiff, and jpg files
def savefigure(buffer):
        plt.figure(figsize=(8, 4))
        oxygenfugacity, H2S_to_total_S, SO2_to_total_S = fluid_s_speciation_curve(buffer)
        plt.plot(oxygenfugacity, H2S_to_total_S, label='H2S', linestyle="-", linewidth=2, color='green')
        plt.plot(oxygenfugacity, SO2_to_total_S, label='SO2', linestyle="-", linewidth=2, color='red')
        plt.title('Sulfur Speciation at ' + str(p_bar) + ' bar and ' + str(t_c) + '°C')
        plt.xlabel('fO2(Δ' + buffer + ')')
        plt.ylabel('Fraction of Total Sulfur')
        plt.legend()
        plt.xlim([-4, 4])
        plt.ylim([0, 1])
        plt.savefig(currentfolder+'/SulfurSpeciation_' + str(t_c) + 'C_' + buffer + '.png')
        plt.savefig(currentfolder + '/SulfurSpeciation_' + str(t_c) + 'C_' + buffer + '.tiff')
        plt.savefig(currentfolder + '/SulfurSpeciation_' + str(t_c) + 'C_' + buffer + '.jpg')


savefigure('NNO')
savefigure('FMQ')


# A function to save data as csv and txt files
def savedata(buffer):
        oxygenfugacity, H2S_to_total_S, SO2_to_total_S = fluid_s_speciation_curve(buffer)
        datatosave = zip(oxygenfugacity, H2S_to_total_S, SO2_to_total_S)
        datatosave = pd.DataFrame(datatosave)
        datatosave.to_csv(currentfolder+'/SulfurSpeciation_' + str(t_c) + 'C_' + buffer + '.csv', sep='\t', header=[('Δ'+buffer),'H2S','SO2'], index=False)
        datatosave.to_csv(currentfolder+'/SulfurSpeciation_' + str(t_c) + 'C_' + buffer + '.txt', sep='\t', header=[('Δ' + buffer), 'H2S', 'SO2'], index=False)


savedata('NNO')
savedata('FMQ')