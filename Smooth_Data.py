#!/usr/bin/python
import numpy
import sys
import matplotlib.pyplot as plt

def parse_TSV(file_handle):
    Header = file_handle.readline().strip().replace("\"","").split()
    Genes = Header[1::]
    data_matrix = {}
    for gene in Genes:
        data_matrix[gene]={}
        data_matrix[gene][0] = [[],[]]
    experiment_counter = 0
    for line in file_handle:
        la = line.strip().split()
        if len(la) == 0:
            experiment_counter+=1
            for gene in Genes:
                data_matrix[gene][experiment_counter] = [[],[]]
        else:
            la = map(float,la)
            time = la[0]
            Values = la[1::]
            for index,value in enumerate(Values):
                data_matrix[Genes[index]][experiment_counter][0].append(time)
                data_matrix[Genes[index]][experiment_counter][1].append(value)
    return data_matrix

perturbation_raw = parse_TSV(open(sys.argv[1]))
perturbation_smoothed = parse_TSV(open(sys.argv[2]))

for gene in perturbation_raw.keys():
    for experiment in perturbation_raw[gene].keys():
        time = perturbation_raw[gene][experiment][0]
        xp = numpy.linspace(time[0],time[-1], 100)
        values = perturbation_raw[gene][experiment][1]
        polynomial = numpy.poly1d(numpy.polyfit(time,values,19))
        smoothed_time = perturbation_smoothed[gene][experiment][0]
        smoothed_values = perturbation_smoothed[gene][experiment][1]
        plt.plot(time,values,'.',smoothed_time,smoothed_values,'.',xp,polynomial(xp),'-')
        plt.title(gene+"_"+str(experiment))
        plt.show()