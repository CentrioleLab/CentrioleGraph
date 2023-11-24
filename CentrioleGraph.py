# Version 20 Apr 2020
# Change from 17 Apr 2020 :
#       Distance between MT wall 0 200 nm
# By Maeva Le Guennec
# To change again centriole width : modify line 34 and 35 :
# Values 300 and 500 according to the centriole width (500 - 300 = 200 nm)


# This is how we can acces the ImageJ API:
# https://imagej.nih.gov/ij/developer/api/allclasses-noframe.html
from ij import IJ, WindowManager
from ij.gui import GenericDialog, Roi, Plot, WaitForUserDialog
from ij.process import ImageConverter
from ij.plugin.frame import RoiManager
from java.awt.event   import ActionListener
from fiji.util.gui import GenericDialogPlus
from ij.measure import ResultsTable
from ij.gui import NonBlockingGenericDialog
from ij.io import OpenDialog
import os

import math, re


lab_prot = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']


class Protein:
    def  __init__(self, start, end, sd_st = 0, sd_e = 0, xpos = 0, \
width = 35, col = '#000000'):
        self.start = start
        self.end = end
        self.sd_st = sd_st
        self.sd_e = sd_e
        self.col = col
        self.width = str(abs(int(width))) # no negative value
        self.xpos = [str((300 - float(xpos)) - float(self.width)/2 ), \
str((500 + float(xpos)) - float(self.width)/2)] 
    def attr_rect(self, plot_max = 900):
        y0 = str(plot_max - self.end)
        height = str(self.end - self.start)
        return {'y0' : y0, 'height' : height, 'x' : self.xpos,\
 'col' : self.col, 'width' : self.width}
    def attr_sd_st(self, plot_max = 900):
        y1 = str((plot_max - self.end) + (self.end - self.start)) # y0 + height
        y2 = str(float(y1) + self.sd_st)
        x = [str(float(self.xpos[0]) + float(self.width)/2), \
str(float(self.xpos[1]) + float(self.width)/2)]
        return {'y1' : y1, 'y2' : y2, 'x' : x, 'col' : self.col}
    def attr_sd_e(self, plot_max = 900):
        y1 = str(plot_max - self.end)
        y2 = str((plot_max - self.end) - self.sd_e) # y0 - sd
        x = [str(float(self.xpos[0]) + float(self.width)/2), \
str(float(self.xpos[1]) + float(self.width)/2)]
        return {'y1' : y1, 'y2' : y2, 'x' : x, 'col' : self.col}


def sdCalc(data):
#    Code from :
#https://www.w3resource.com/python-exercises/math/python-math-exercise-57.php
    n = len(data)
    if n <= 1:
        return 0.0
    mean, sd = avgCalc(data), 0.0
    # calculate stan. dev.
    for el in data:
        sd += (float(el) - mean)**2
    sd = math.sqrt(sd / float(n-1))
    return round(sd, 2)


def avgCalc(ls):
#    Code from :
#https://www.w3resource.com/python-exercises/math/python-math-exercise-57.php
    n, mean = len(ls), 0.0
    if n <= 1:
        return ls[0]
    # calculate average
    for el in ls:
        mean = mean + float(el)
    mean = mean / float(n)
    return round(mean, 2)


def writeHeader900(filout):
# Y axis goes to 800
    filout.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
    filout.write('<svg xmlns="http://www.w3.org/2000/svg" \
xmlns:xlink="http://www.w3.org/1999/xlink" width="700px" height="1100px" \
version="1.1">\n')
    filout.write('<rect width="700" height="1100" fill="#ffffff" />\n')
    filout.write('<polyline points="175, 50, 175, 900, 630, 900" \n')
    filout.write('fill="white" stroke="#000000" stroke-width="3" /> \
<!--XY axis, starts from x,y = 175,900--> \n')
    filout.write('<!--Tickmarks along Y axis--> \n')
    filout.write('<line x1="150" y1="700" x2="175" y2="700" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<line x1="150" y1="500" x2="175" y2="500" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<line x1="150" y1="300" x2="175" y2="300" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<line x1="150" y1="100" x2="175" y2="100" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<!--Half-Tickmarks along Y axis--> \n')
    filout.write('<line x1="165" y1="800" x2="175" y2="800" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<line x1="165" y1="600" x2="175" y2="600" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<line x1="165" y1="400" x2="175" y2="400" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<line x1="165" y1="200" x2="175" y2="200" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<!--Y-labels-->\n')
    filout.write('<text x="75" y="710" fill="#000000" font-size="28" \
font-family="\'Raleway\', sans-serif" > 200 </text>\n')
    filout.write('<text x="75" y="510" fill="#000000" font-size="28" \
font-family="\'Raleway\', sans-serif" > 400 </text>\n')
    filout.write('<text x="75" y="310" fill="#000000" font-size="28" \
font-family="\'Raleway\', sans-serif" > 600 </text>\n')
    filout.write('<text x="75" y="110" fill="#000000" font-size="28" \
font-family="\'Raleway\', sans-serif" > 800 </text>\n')
    filout.write('<text x="40" y="410" fill="#000000" font-size="34" \
font-family="\'Raleway\', sans-serif" \n')
    filout.write('    style="writing-mode: tb;" transform="rotate(180,40,470)">\
 Position (nm) </text> \n')



def writeHeader800(filout):
# Y axis goes to 700
    filout.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
    filout.write('<svg xmlns="http://www.w3.org/2000/svg" \
xmlns:xlink="http://www.w3.org/1999/xlink" width="700px" height="1000px"\
 version="1.1">\n')
    filout.write('<rect width="700" height="1000" fill="#ffffff" />\n')
    filout.write('<polyline points="175, 50, 175, 800, 630, 800" \n')
    filout.write('fill="white" stroke="#000000" stroke-width="3" /> \
<!--XY axis, starts from x,y = 175,900--> \n')
    filout.write('<!--Tickmarks along Y axis--> \n')
    filout.write('<line x1="150" y1="600" x2="175" y2="600" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<line x1="150" y1="400" x2="175" y2="400" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<line x1="150" y1="200" x2="175" y2="200" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<!--Half-Tickmarks along Y axis--> \n')
    filout.write('<line x1="165" y1="700" x2="175" y2="700" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<line x1="165" y1="500" x2="175" y2="500" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<line x1="165" y1="300" x2="175" y2="300" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<line x1="165" y1="100" x2="175" y2="100" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<!--Y-labels-->\n')
    filout.write('<text x="75" y="610" fill="#000000" font-size="28" \
font-family="\'Raleway\', sans-serif" > 200 </text>\n')
    filout.write('<text x="75" y="410" fill="#000000" font-size="28" \
font-family="\'Raleway\', sans-serif" > 400 </text>\n')
    filout.write('<text x="75" y="210" fill="#000000" font-size="28" \
font-family="\'Raleway\', sans-serif" > 600 </text>\n')
    filout.write('<text x="40" y="310" fill="#000000" font-size="34" \
font-family="\'Raleway\', sans-serif" \n')
    filout.write('    style="writing-mode: tb;" transform="rotate(180,40,470)">\
 Position (nm) </text> \n')


def writeHeader700(filout):
# Y axis goes to 600
    filout.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
    filout.write('<svg xmlns="http://www.w3.org/2000/svg" \
xmlns:xlink="http://www.w3.org/1999/xlink" width="700px" height="900px" \
version="1.1">\n')
    filout.write('<rect width="700" height="900" fill="#ffffff" />\n')
    filout.write('<polyline points="175, 50, 175, 700, 630, 700" \n')
    filout.write('fill="white" stroke="#000000" stroke-width="3" /> \
<!--XY axis, starts from x,y = 175,900--> \n')
    filout.write('<!--Tickmarks along Y axis--> \n')
    filout.write('<line x1="150" y1="500" x2="175" y2="500" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<line x1="150" y1="300" x2="175" y2="300" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<line x1="150" y1="100" x2="175" y2="100" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<!--Half-Tickmarks along Y axis--> \n')
    filout.write('<line x1="165" y1="600" x2="175" y2="600" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<line x1="165" y1="400" x2="175" y2="400" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<line x1="165" y1="200" x2="175" y2="200" stroke="#000000" \
stroke-width="3"/>\n')
    filout.write('<!--Y-labels-->\n')
    filout.write('<text x="75" y="510" fill="#000000" font-size="28" \
font-family="\'Raleway\', sans-serif" > 200 </text>\n')
    filout.write('<text x="75" y="310" fill="#000000" font-size="28" \
font-family="\'Raleway\', sans-serif" > 400 </text>\n')
    filout.write('<text x="75" y="110" fill="#000000" font-size="28" \
font-family="\'Raleway\', sans-serif" > 600 </text>\n')
    filout.write('<text x="40" y="410" fill="#000000" font-size="34" \
font-family="\'Raleway\', sans-serif" \n')
    filout.write('    style="writing-mode: tb;" transform="rotate(180,40,470)">\
 Position (nm) </text> \n')


def writeProtein(filout, protein, y_max, width = '30'):
        p_rect = protein.attr_rect(int(y_max))
        p_sd_st = protein.attr_sd_st(int(y_max))
        p_sd_e = protein.attr_sd_e(int(y_max))
        dist = '0'
        sd_col = protein.col # colored sd bars
        filout.write('    <rect width="' + str(p_rect['width']) + \
'" height=" '+ p_rect['height'] + '" fill="' + protein.col +'" x="' + \
str(p_rect['x'][0]) + '" y="' + p_rect['y0'] + '"/>\n')
        filout.write('z<rect width="' + str(p_rect['width']) + \
'" height=" '+ p_rect['height'] + '" fill="' + protein.col +'" x="' + \
str(p_rect['x'][1]) + '" y="' + p_rect['y0'] + '"/>\n')
# SD bars (bottom, starting point)
        filout.write('    <line x1="' + p_sd_st['x'][0] + '" y1="' + \
p_sd_st['y1'] + '" x2="' + p_sd_st['x'][0] + '" y2="' + p_sd_st['y2'] +\
'" stroke="' + sd_col + '" stroke-width="3" stroke-dasharray="10,' + dist +'"/>\n')
        filout.write('    <line x1="' + p_sd_st['x'][1] + '" y1="' + \
p_sd_st['y1'] + '" x2="' + p_sd_st['x'][1] + '" y2="' + p_sd_st['y2'] + \
'" stroke="' + sd_col + '" stroke-width="3" stroke-dasharray="10,' + dist +'"/>\n')
# SD bars (top, end point)
        filout.write('    <line x1="' + p_sd_e['x'][0] + '" y1="' + \
p_sd_e['y1'] + '" x2="' + p_sd_e['x'][0] + '" y2="' + p_sd_e['y2'] + \
'" stroke="' + sd_col + '" stroke-width="3" stroke-dasharray="10,' + dist +'"/>\n')
        filout.write('    <line x1="' + p_sd_e['x'][1] + '" y1="' + \
p_sd_e['y1'] + '" x2="' + p_sd_e['x'][1] + '" y2="' + p_sd_e['y2'] + \
'" stroke="' + sd_col + '" stroke-width="3" stroke-dasharray="10,' + dist +'"/>\n')
# SD bars horizontal bottom/start
        filout.write('    <line x1="' + protein.xpos[0] + '" y1="' + \
p_sd_st['y2'] + '" x2="' + str(float(protein.xpos[0]) + int(protein.width)) + \
'" y2="' + p_sd_st['y2'] + \
'" stroke="' + sd_col + '" stroke-width="3" stroke-dasharray="10,' + dist +'"/>\n')
        filout.write('    <line x1="' + protein.xpos[1] + '" y1="' + \
p_sd_st['y2'] + '" x2="' + str(float(protein.xpos[1]) + int(protein.width)) + \
'" y2="' + p_sd_st['y2'] + \
'" stroke="' + sd_col + '" stroke-width="3" stroke-dasharray="10,' + dist +'"/>\n')
# SD bars horizontal top/end
        filout.write('    <line x1="' + protein.xpos[0] + '" y1="' + \
p_sd_e['y2'] + '" x2="' + str(float(protein.xpos[0]) + int(protein.width)) + \
'" y2="' + p_sd_e['y2'] + \
'" stroke="' + sd_col + '" stroke-width="3" stroke-dasharray="10,' + dist +'"/>\n')
        filout.write('    <line x1="' + protein.xpos[1] + '" y1="' + \
p_sd_e['y2'] + '" x2="' + str(float(protein.xpos[1]) + int(protein.width)) + \
'" y2="' + p_sd_e['y2'] + \
'" stroke="' + sd_col + '" stroke-width="3" stroke-dasharray="10,' + dist +'"/>\n')


def writeFile(filename, list_prots, plotname = 'Plot Name', height = 900):
    filout = open(filename, 'w')
    if height == '900' :
        writeHeader900(filout)
        filout.write('    <text x="400" y="1000" fill="#000000" \
font-size="34" font-family="\'Raleway\', sans-serif" text-anchor="middle" >')
    elif height == '800' :
        writeHeader800(filout)
        filout.write('    <text x="400" y="900" fill="#000000" \
font-size="34" font-family="\'Raleway\', sans-serif" text-anchor="middle" >')
    else:
        writeHeader700(filout)
        filout.write('    <text x="400" y="800" fill="#000000" \
font-size="34" font-family="\'Raleway\', sans-serif" text-anchor="middle" >')
    filout.write(plotname + '</text> \n')
# Rectangles
    for p in list_prots:
        writeProtein(filout, p, height)
    filout.write('</svg>')
    filout.close()

def readCSVFile(filename):
    filin = open(filename, 'r')
    header = filin.readline()
    if len(header.split(','))%2 == 0 :
        print('Error, wrong number of columns !')
        print('You should have "Label" followed by 2 columns for' +\
' each protein ("start" and "end")  ')
        return False
    elif len(header.split(',')) > 10 :
        print('Too many columns')
        print('There are more than ten proteins in your file.' )
        return False
    nb_p =  int((len(header.split(',')) - 1)/2 )
    dico_lists_st = {} 
    dico_lists_en = {}
    for k in range(nb_p):
        dico_lists_st[lab_prot[k]] = []
        dico_lists_en[lab_prot[k]] = []
    for line in filin:
        split_line = line.strip().split(',')
        i = 0
        for j in range(1, nb_p + 2, 2):
            dico_lists_st[lab_prot[i]].append(\
round(float(split_line[j]), 2))
            dico_lists_en[lab_prot[i]].append(\
round(float(split_line[j + 1]), 2))
            i += 1
    return ({'start_data' : dico_lists_st,'end_data' : dico_lists_en})


def shiftValues(dico_data, ref):
    dico_shifted = {}
    dico_shifted['start_data'] = {}
    dico_shifted['end_data'] = {}
    for prot in dico_data['start_data'].keys() :
        dico_shifted['start_data'][prot] = []
        dico_shifted['end_data'][prot] = []
    for prot in dico_data['start_data'].keys() : # for each protein
        for i in range(len(dico_data['start_data'][prot])): # for each value
            val_st = dico_data['start_data'][prot][i] - \
dico_data['start_data'][ref][i]
            val_en = dico_data['end_data'][prot][i] - \
dico_data['start_data'][ref][i]
            dico_shifted['start_data'][prot].append(val_st)
            dico_shifted['end_data'][prot].append(val_en)
    return dico_shifted


def rescValues(dico_data, ref):
    dico_resc = {}
    dico_resc['start_data'] = {}
    dico_resc['end_data'] = {}
    for prot in dico_data['start_data'].keys() :
        dico_resc['start_data'][prot] = []
        dico_resc['end_data'][prot] = []
    mean_val = avgCalc(dico_data['end_data'][ref]) # avg length value
    for prot in dico_data['start_data'].keys() :
        for i in range(len(dico_data['start_data'][prot])): # for each value
            if prot != ref:
                new_st = dico_data['start_data'][prot][i] * \
(mean_val/dico_data['end_data'][ref][i])
                new_en = dico_data['end_data'][prot][i] * \
(mean_val/dico_data['end_data'][ref][i])
                dico_resc['start_data'][prot].append(new_st)
                dico_resc['end_data'][prot].append(new_en)
            else :
                dico_resc['start_data'][prot].append(\
dico_data['start_data'][prot][i])
                dico_resc['end_data'][prot].append(\
dico_data['end_data'][prot][i])
    return dico_resc


def statsList(dico_lists):
    dico_stats = {}
    for i in dico_lists.keys() :
        dico_stats[i] = {'mean' : avgCalc(dico_lists[i]), \
'sd' : sdCalc(dico_lists[i])}
    return dico_stats


def createProteins(dico_stats, list_xpos, list_wid, list_col, ref):
    list_proteins = []
    nb_prot = len(dico_stats['start_stats'])
    i = 0
    for label in dico_stats['start_stats'].keys():
        protein = Protein(dico_stats['start_stats'][label]['mean'],\
dico_stats['end_stats'][label]['mean'], dico_stats['start_stats'][label]['sd'],\
dico_stats['end_stats'][label]['sd'], list_xpos[i], list_wid[i], list_col[i])
        list_proteins.append(protein)
        i += 1
    return(list_proteins)


def txtToList(names, nb_prot):
    names = names.replace(' ', '')
    names = names.replace('\t', '')
    list_val = names.split(',')
    while len(list_val) < nb_prot:  # fill the colnames list until enough values
            list_val.append(list_val[len(list_val) - 1])
    return list_val

def colOK(col_text):
    col_text = col_text.replace(' ', '')
    col_text = col_text.replace('\t', '')
    list_col = col_text.split(',')
    for i in list_col :# check if i is an hexadecimal code
        match = re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', i) 
        if match:                      
            continue
        else:
            return False
    return True

def numericOK(text):
    text = text.replace(' ', '')
    text = text.replace('\t', '')
    pos = text.split(',')
    for i in pos : # check if i is an float
        match = re.search(r'^[-+]{0,1}[0-9]+(\.[0-9]+){0,1}$', i)
        if match:                      
            continue
        else:
            return False
    return True


def pixToNm(dico, factor):
    dico_nm = {}
    for key in dico.keys():
        dico_nm[key] = dico[key]
        for key2 in dico[key].keys():
            dico_nm[key][key2] = dico[key][key2]
            for key3 in dico[key][key2].keys():
                dico_nm[key][key2][key3] = dico[key][key2][key3]*factor
    return dico_nm


def testEntriesOK(pos_text, width_text, col_text):
    if (numericOK(pos_text)) :
        if (numericOK(width_text)) :
            if colOK(col_text) :
                return True
            else :
                return 'col'
        else :
            return 'wid'
    else :
        return 'pos'


def mainPipeline_resc(filin_name, filout_name, plot_name,\
factor, ref, pos, wid, col, maxY):
    dico_columns = readCSVFile(filin_name)
    dico_shifted = shiftValues(dico_columns, ref)
    dico_resc = rescValues(dico_shifted, ref)
    dico_stats = {'start_stats':statsList(dico_resc['start_data']),\
'end_stats':statsList(dico_resc['end_data'])}
    dico_stats_nm = pixToNm(dico_stats, factor)
    nb_prots = len(dico_columns.keys())
    list_col = txtToList(col, nb_prots)
    list_pos = txtToList(pos, nb_prots)
    list_wid = txtToList(wid, nb_prots)
    list_proteins = createProteins(dico_stats_nm, \
list_pos, list_wid, list_col, ref)
    writeFile(filout_name, list_proteins, plot_name, maxY)


#def mainPipeline_noresc(filin_name, filout_name, plot_name,\
#factor, ref, pos, col, maxY):
#    dico_columns = readCSVFile(filin_name)
#    dico_shifted = shiftValues(dico_columns, ref)
#    dico_stats = {'start_stats':statsList(dico_shifted['start_data']),\
#'end_stats':statsList(dico_shifted['end_data'])}
#    dico_stats_nm = pixToNm(dico_stats, factor)
#    list_proteins = createProteins(dico_stats_nm, list_pos, list_col, ref)
#    writeFile(filout_name, list_proteins, plot_name, maxY)



def menuWindow(filin, filout, fact, ref, pos, width, col, title, maxY):
#    win = NonBlockingGenericDialog("Options Choices")
    win = GenericDialogPlus("Position plot - v.17Apr2020")
# Output file
    win.addMessage('---------------------------------------------------------  \
Data info --------------------------------------------------------------------')
    win.addFileField(\
"Enter the input file (.csv, the first column must be a ''Label'' column):",\
filin, 40)
    win.addFileField("Enter the output file (a .svg file):",\
filout, 40)
# Proteins info 
    win.addNumericField('Factor to apply to convert values to nm : ', 1.000, 3)
    win.addMessage('Indicate which protein should be used as for \
length reference.')
    win.addChoice('Ex: If tubulin is the first protein, select 1', \
['1', '2','3','4','5','6','7','8','9','10'], ref)
    win.addMessage('---------------------------------------------------------  \
Proteins appearance ----------------------------------------------------------')
    win.addMessage('For the following parameters, you can use only 1 value \
each that will be common for all proteins.')
    win.addMessage('If the parameters are different for each protein, values \
entered must be separated by ",".')
    win.addMessage('Distance of between the protein and\
 the microtubule wall in nm :')
    win.addStringField('(< 0 : internal protein, > 0 : external protein)',\
pos, 20)
    win.addStringField('Thickness of the protein bar (in nm) : ',\
width, 20)
    win.addStringField('Colors as hexadecimal code (#000000) : ', \
col, 20)
# Plot features 
    win.addMessage('---------------------------------------------------------  \
Plot features ----------------------------------------------------------------')
    win.addStringField("Plot title:", title, 20)
    win.addChoice('Select the Y axis maximum:', ['600', '700', '800'], maxY)
# Show dialog and save value 
    win.showDialog()
    input_name = win.getNextString()
    output_name = win.getNextString()
    factor = win.getNextNumber() 
    prot_ref = str(win.getNextChoiceIndex() + 1) # index goes from 0 to 9
    prot_pos = win.getNextString()
    prot_wid = win.getNextString()
    prot_col = win.getNextString()
    plot_name = win.getNextString()
    Ymax_index = win.getNextChoiceIndex()
    new_maxY = ['600', '700', '800'][Ymax_index] # for menu
    y_axis = ['700', '800', '900'][Ymax_index] # for mainpipeline
    if win.wasCanceled():
        return False
    elif win.wasOKed():
        if (testEntriesOK(prot_pos, prot_wid, prot_col) == True) :
            if (os.path.exists(output_name)): # Check output exists or not
                war = GenericDialog("Warning")
                war.addMessage('File Exists and will be replaced.')
                war.enableYesNoCancel('Change name', 'Keep this name')
                war.hideCancelButton()
                war.showDialog()
                if war.wasOKed():
                    menuWindow(input_name, output_name, factor,\
prot_ref, prot_pos, prot_wid, prot_col, plot_name, new_maxY)
                else :
                    mainPipeline_resc(input_name, output_name, plot_name,\
factor, prot_ref, prot_pos, prot_wid, prot_col, y_axis)
            else :
                mainPipeline_resc(input_name, output_name, plot_name,\
factor, prot_ref, prot_pos, prot_wid, prot_col, y_axis)
            menuWindow(input_name, output_name, factor, prot_ref, \
prot_pos, prot_wid, prot_col, plot_name, new_maxY)
        elif (testEntriesOK(prot_pos, prot_wid, prot_col) == 'wid') :
            wait = WaitForUserDialog("Error in width parameters.")
            wait.show()
            menuWindow(input_name, output_name, factor, prot_ref, \
prot_pos, prot_wid, prot_col, plot_name, new_maxY)
        elif (testEntriesOK(prot_pos, prot_wid, prot_col) == 'col') :
            wait = WaitForUserDialog("Error in colors parameters.")
            wait.show()
            menuWindow(input_name, output_name, factor, prot_ref, \
prot_pos, prot_wid, prot_col, plot_name, new_maxY)
        elif (testEntriesOK(prot_pos, prot_wid, prot_col) == 'pos') :
            wait = WaitForUserDialog("Error in positions parameters.")
            wait.show()
            menuWindow(input_name, output_name, factor, prot_ref, \
prot_pos, prot_wid, prot_col, plot_name, new_maxY)




#def inputFile():
#    srcDir = None
#    while srcDir == None:
#        a = OpenDialog('Select an input file')
#        srcDir = a.getDirectory()
#    # User canceled the dialog
#    return (a.getFileName(), srcDir)

if __name__ in ['__builtin__','__main__']:
    folder = IJ.getDirectory("current")
    results = menuWindow(folder, folder + 'output.svg',\
1, '1', '0, -15', '35, 35',\
"#808080, #ff0080", 'Title', '600')







