#!/usr/bin/env python3
#
from __future__ import print_function
import textwrap
import datetime
import sys
from decimal import *  # noqa

# inputs: *.pen file

# outputs: *.mc file for MCNP5 input, *.ref file

# description of *.ref file
# -------------------------

#   The reference file is similar to that of a PENTRAN .flx file,
#   with all information except for the flux,
#   which is obtained in another post-processing script.
#   It has the format Group Coarse Dom.Matl x y z Flux(Blank),
#   which makes the MCNP file appear in a format that is similar
#   to that of a PENTRAN file

# handling of eigenvalue problems
# -------------------------------

#   If a 'kcode' or 'kcode.txt' is found in the same directory
#   The file's contents are spliced into the MCNP5 data section


#class hashabledict(dict):
#    def __hash__(self):
#        return hash(tuple(sorted(self.items())))

### defining an instance class

class Instance:
    def __init__(self, value):
        self.value = value

    def stringArgs(self, prestring):
        """Method docstring."""
        return prestring+str(self.value)

instances = {}


def generate_instance(name, value):
    myName = name.strip()
    myValue = value
    instances[myName] = Instance(myValue)

### end of instance class

# auxilary functions


def pretty_print(var_str):
    print(var_str+':', len(globals()[var_str]), globals()[var_str])
    print()


def table_print(dictionary, var_str):
    table = []
    name_list = []
    for item in var_str.split():
        name_list.append(item)
        table.append(globals()[dictionary][item])
        # table_length=len(globals()[dictionary][item])
    print_table = zip(*table)

    print(' C.M.', end="")
    for item in name_list:
        print(item.rjust(11), end="")
    print()

    fm_count = 0
    for item in print_table:
        fm_count += 1
        print('{0:5d}'.format(fm_count), end="")
        for sub_item in item:
            #print(subitem)
            print('   {0:8.{1}f}'.format(float(sub_item), 4), end="")
        print()
    print()


def fido_grab(line, fido_per_cm):
    """ fido_per_cm needs to be instantiated before calling """
    count = 0
    if "nmattp" in line:

        line_list = line.split('=')  # clean out the '=' sign
        line_list.pop(0)           # pop 'nmattp' 1st item off

        string = line_list.pop()
        line_list = string.split()   # the number string should be split
        # cmCount = int(line_list.pop(0))
        # lump the fine mesh fido values back together in one string
        fido_string = ' '.join(line_list)
        fido_per_cm.append(fido_string)

        count += 1

    return fido_per_cm


def fine_check(line):
    """ generate instances for variables having the word 'fine' """
    if "fine" in line:
        line_list = line.split('=')
        val = line_list.pop().replace("\n", "")
        type = line_list.pop()
        generate_instance(type, val)


def mesh_boundary_check(line):
    """ generate instances for variables having the word 'mesh' """
    if "mesh" in line:
        line_list = line.split('=')
        val = line_list.pop().replace("\n", "")
        type = line_list.pop()
        generate_instance(type, val)


def cm_check(f, line):
    if "maxgcm" in line:
        values = next(f)               # grab the next line !
        values = values.split()
        maxgcm = int(values[2])        # pull the 3rd value for maxgcm
        generate_instance("maxgcm", maxgcm)


def grp_matl_check(f, line):
    if "maxmat" in line:
        values = next(f)               # grab the next line !
        values = values.split()
        maxgrp = int(values[0])
        maxmat = int(values[4])        # pull the 5th value for maxgcm
        generate_instance("maxmat", maxmat)
        generate_instance("maxgrp", maxgrp)


def bdy_check(f, line):
    if "ibback" in line:
        values = line.split('=')
        ibback = values[1]
        generate_instance("ibback", ibback)
    if "ibfrnt" in line:
        values = line.split('=')
        ibfrnt = values[1]
        generate_instance("ibfrnt", ibfrnt)
    if "jbeast" in line:
        values = line.split('=')
        jbeast = values[1]
        generate_instance("jbeast", jbeast)
    if "jbwest" in line:
        values = line.split('=')
        jbwest = values[1]
        generate_instance("jbwest", jbwest)
    if "kbsout" in line:
        values = line.split('=')
        kbsout = values[1]
        generate_instance("kbsout", kbsout)
    if "kbnort" in line:
        values = line.split('=')
        kbnort = values[1]
        generate_instance("kbnort", kbnort)


def calc_fm_per_cm():
    """ Calculate Fine Mesh per Coarse Mesh """
    print("\nCalculating number of fine meshes per coarse mesh:")
    print("-"*80)
    # maxgcm = instances["maxgcm"].value
    ixL = instances["ixfine"].value.split()
    jyL = instances["jyfine"].value.split()
    kzL = instances["kzfine"].value.split()
    ixL = [int(val) for val in ixL]
    jyL = [int(val) for val in jyL]
    kzL = [int(val) for val in kzL]

    # the next one is a trick!
    # result_string = ' '.join(
    # [str(index + 1) + "-->" + str(i*j*k)
    #     for index, (i, j, k) in enumerate(list(zip(ixL, jyL, kzL)))])

    result_list = \
        [str(index + 1) + "-->" + str(i*j*k)
            for index, (i, j, k) in enumerate(list(zip(ixL, jyL, kzL)))]

    ### just enforcing a prettier output
    formatted_string = ""
    for item in result_list:
            formatted_string += str(item.ljust(14))
    print(textwrap.fill(formatted_string, width=84)+"\n")

    print("Total FM in all CM-->" +
          str(sum([i*j*k for i, j, k in zip(ixL, jyL, kzL)])))
    fm_per_cm = sum([i*j*k for i, j, k in zip(ixL, jyL, kzL)])

    fm_attributes = dict(fm_per_cm=fm_per_cm,
                         xFm=ixL,
                         yFm=jyL,
                         zFm=kzL)
    return fm_attributes


def calc_boundary_per_cm(xFm, yFm, zFm):
    """ Collect physical boundaries into six 1D arrays for each coarse mesh """
    print("\nCalculating physical boundaries per coarse mesh:")
    print("-"*80)
    xMinCm = []
    xMaxCm = []
    xFmDimL = []
    yMinCm = []
    yMaxCm = []
    yFmDimL = []
    zMinCm = []
    zMaxCm = []
    zFmDimL = []
    xmeshL = instances["xmesh"].value.split()
    ymeshL = instances["ymesh"].value.split()
    zmeshL = instances["zmesh"].value.split()
    x_dim_span = len(xmeshL) - 1
    y_dim_span = len(ymeshL) - 1
    z_dim_span = len(zmeshL) - 1
    maxgcm = x_dim_span * y_dim_span * z_dim_span
    print("Independent calculation of maxgcm-->"+str(maxgcm))
    idx = -1
    for k in range(z_dim_span):
        for j in range(y_dim_span):
            for i in range(x_dim_span):
                idx += 1
                xMinCm.append(float(xmeshL[i]))
                xMaxCm.append(float(xmeshL[i+1]))
                yMinCm.append(float(ymeshL[j]))
                yMaxCm.append(float(ymeshL[j+1]))
                zMinCm.append(float(zmeshL[k]))
                zMaxCm.append(float(zmeshL[k+1]))
                xFmDimL.append(Decimal((xMaxCm[idx] - xMinCm[idx])) /
                               Decimal(xFm[idx]))
                yFmDimL.append(Decimal((yMaxCm[idx] - yMinCm[idx])) /
                               Decimal(yFm[idx]))
                zFmDimL.append(Decimal((zMaxCm[idx] - zMinCm[idx])) /
                               Decimal(zFm[idx]))

    ixfL = instances["ixfine"].value.split()
    jyfL = instances["jyfine"].value.split()
    kzfL = instances["kzfine"].value.split()

    ixfL = [int(i) for i in ixfL]
    jyfL = [int(i) for i in jyfL]
    kzfL = [int(i) for i in kzfL]

    fm_per_cm = [p * q * r for p, q, r in zip(ixfL, jyfL, kzfL)]

    cm_attributes = dict(
        xMinCm=xMinCm, xMaxCm=xMaxCm,
        yMinCm=yMinCm, yMaxCm=yMaxCm,
        zMinCm=zMinCm, zMaxCm=zMaxCm,
        xFmDimL=xFmDimL, yFmDimL=yFmDimL, zFmDimL=zFmDimL,
        ixfL=ixfL, jyfL=jyfL, kzfL=kzfL,
        fm_per_cm=fm_per_cm)

    return cm_attributes


def produce_surface_rpp(cm_attributes):

    xMinCm = cm_attributes['xMinCm']
    xMaxCm = cm_attributes['xMaxCm']
    yMinCm = cm_attributes['yMinCm']
    yMaxCm = cm_attributes['yMaxCm']
    zMinCm = cm_attributes['zMinCm']
    zMaxCm = cm_attributes['zMaxCm']
    xFmDimL = cm_attributes['xFmDimL']
    yFmDimL = cm_attributes['yFmDimL']
    zFmDimL = cm_attributes['zFmDimL']

    rpp_list = []
    for idx in range(0, len(xMinCm)):

        rpp_desc_str = "c " + " "*14 + "FM dimensions for CM:" + \
            str(int(idx + 1))+"\n"

        rpp_string = rpp_desc_str + str(int(idx + 1)).ljust(5) + " rpp "

        str_min = '%.15E' % xMinCm[idx]
        str_max = '%.15E' % (Decimal(xMinCm[idx]) + xFmDimL[idx])
        rpp_string += str_min.ljust(20) + "\n" + " "*10 + str_max.ljust(20) +\
            " $ x-min x-max\n" + " "*10

        str_min = '%.15E' % yMinCm[idx]
        str_max = '%.15E' % (Decimal(yMinCm[idx]) + yFmDimL[idx])
        rpp_string += str_min.ljust(20) + "\n" + " "*10 + str_max.ljust(20) + \
            " $ y-min y-max"

        rpp_list.append(rpp_string)

        str_min = '%.15E' % zMinCm[idx]
        str_max = '%.15E' % (Decimal(zMinCm[idx]) + zFmDimL[idx])
        rpp_string = " "*10 + str_min.ljust(20) + "\n" + " "*10 + \
            str_max.ljust(20) + " "

        rpp_string += "$ z-min z-max\nc"

        rpp_list.append(rpp_string)

    for idx in range(0, len(xMinCm)):

        rpp_desc_str = "c " + " "*14 + "CM dimensions for CM:" + \
            str(int(idx + 1)) + "\n"

        rpp_string = rpp_desc_str + str(int(idx + 1 + 10000)).ljust(5) + \
            " rpp "

        str_min = '%.14e' % xMinCm[idx]
        str_max = '%.14e' % xMaxCm[idx]
        rpp_string += str_min.ljust(20) + "\n" + " "*10 + \
            str_max.ljust(20) + " $ x-min x-max\n"+" " * 10

        str_min = '%.14e' % yMinCm[idx]
        str_max = '%.14e' % yMaxCm[idx]
        rpp_string += str_min.ljust(20) + "\n" + " "*10 + str_max.ljust(20) + \
            " $ y-min y-max"

        rpp_list.append(rpp_string)

        str_min = '%.14e' % zMinCm[idx]
        str_max = '%.14e' % zMaxCm[idx]
        rpp_string = " "*10 + str_min.ljust(20) + "\n" + " "*10 + \
            str_max.ljust(20) + " "

        rpp_string += "$ z-min z-max\nc"

        rpp_list.append(rpp_string)

    # the global rpp container
    xmin = min(xMinCm)
    xmax = max(xMaxCm)
    ymin = min(yMinCm)
    ymax = max(yMaxCm)
    zmin = min(zMinCm)
    zmax = max(zMaxCm)

    xmin_rpp_string = ""
    xmax_rpp_string = ""
    ymin_rpp_string = ""
    ymax_rpp_string = ""
    zmin_rpp_string = ""
    zmax_rpp_string = ""

    if "1" in instances["ibback"].value:
        xmin_rpp_string = "*"
    if "1" in instances["ibfrnt"].value:
        xmax_rpp_string = "*"
    if "1" in instances["jbeast"].value:
        ymin_rpp_string = "*"
    if "1" in instances["jbwest"].value:
        ymax_rpp_string = "*"
    if "1" in instances["kbsout"].value:
        zmin_rpp_string = "*"
    if "1" in instances["kbnort"].value:
        zmax_rpp_string = "*"

    rpp_list.append("c asterisk indicates reflected boundaries")
    xmin_rpp_string += "99991 px " + "{0:.8f}".format(xmin - 0.00000001)
    rpp_list.append(xmin_rpp_string)
    xmax_rpp_string += "99992 px " + "{0:.8f}".format(xmax + 0.00000001)
    rpp_list.append(xmax_rpp_string)
    ymin_rpp_string += "99993 py " + "{0:.8f}".format(ymin - 0.00000001)
    rpp_list.append(ymin_rpp_string)
    ymax_rpp_string += "99994 py " + "{0:.8f}".format(ymax + 0.00000001)
    rpp_list.append(ymax_rpp_string)
    zmin_rpp_string += "99995 pz " + "{0:.8f}".format(zmin - 0.00000001)
    rpp_list.append(zmin_rpp_string)
    zmax_rpp_string += "99996 pz " + "{0:.8f}".format(zmax + 0.00000001)
    rpp_list.append(zmax_rpp_string)
    void_rpp_string = "99990 rpp " + str(xmin) + " " + str(xmax) + " " + \
        str(ymin) + " " + str(ymax) + " " + \
        str(zmin) + " " + str(zmax)
    rpp_list.append("c global rpp")
    rpp_list.append(void_rpp_string)

    return rpp_list

# expandFido routines


def expand_repeats(fido_string):
    expand_string = ""
    rCount = 0
    fido_list = fido_string.split()
    for item in fido_list:
        if "R" in item:
            rCount += 1
            repeats, matID = item.split("R")
            for i in range(0, int(repeats)):
                expand_string += matID + " "
        else:
            expand_string += item + " "

    expand_repeat_string = expand_string
    return expand_repeat_string, rCount


def expand_q_mults(fido_string):
    expand_string = ""
    q_count = 0
    fido_list = fido_string.split()

    for index in range(len(fido_list)):
        if "Q" in fido_list[index] and q_count == 0:
            q_count += 1
            mults, seq_length = fido_list[index].split("Q")
            for i in range(0, int(mults)):
                for j in range(0, int(seq_length)):
                    expand_string += fido_list[index - int(seq_length) + j] \
                        + " "
        else:
            expand_string += fido_list[index] + " "

    expand_mults_string = expand_string
    return expand_mults_string, q_count


def mcnp_r_fido(full_string):
    full_list = full_string.split()

    count = 0
    gcount = 0
    mcnp_string = full_list[0] + " "

    for index in range(1, len(full_list)):
        if int(full_list[index]) == int(full_list[index - 1]):
            count += 1
        else:
            if count >= 1:
                mcnp_string += str(count) + "R " + str(full_list[index])+" "
                gcount += count
            else:
                mcnp_string += full_list[index] + " "
            gcount += 1
            count = 0
    if count >= 1:
                mcnp_string += str(count) + "R " + " "
                gcount += count + 1

    mcnp_string = \
        textwrap.fill(mcnp_string, width=80,
                      initial_indent='     ',
                      subsequent_indent='     ')
    return mcnp_string


def process_pfm_string(fido_string):
    # expand all R sequences
    expand_repeat_string, rCount = expand_repeats(fido_string)

    # Q sequence repeater has to be iterated after each Q expansion
    # so that no 'Q' items are left in the input
    q_count = 1
    while q_count == 1:
        expand_mults_string, q_count = expand_q_mults(expand_repeat_string)
        expand_repeat_string = expand_mults_string

    mcnp_string = mcnp_r_fido(expand_repeat_string)
    return mcnp_string, expand_repeat_string


# cell card production routines
def produce_matl_cell():
    maxmat = instances["maxmat"].value
    maxmat = int(maxmat)
    matlCell = []
    for i in range(1, maxmat + 1):
        string = str(i).ljust(5)
        string += str(i).ljust(4) + "1.0  -99990   u="
        string += str(i).ljust(4) + "  imp:n=1 $ material " + str(i)
        matlCell.append(string)
    return matlCell


def pretty_print_instances(instances):
    print()
    itemString = ""
    for item in instances:
        prestring = str(item) + "--> "
        if "mesh" in str(item):
            meshString = instances[item].stringArgs(prestring)+" "
            print(textwrap.fill(meshString, width=80, subsequent_indent=' '))
        elif "fine" in str(item):
            fineString = instances[item].stringArgs(prestring)+" "
            print(textwrap.fill(fineString, width=80, subsequent_indent=' '))
        else:
            itemString += instances[item].stringArgs(prestring) + " "

    print(textwrap.fill(itemString, width=80, subsequent_indent=' '))
    print()


def header():
    print(" ")
    print("                          MESH2MCNP")
    print("                         Version 1.0e")
    print("                        By  K. Manalo")
    print("                     Georgia Tech / CRITCEL")
    print("                          Nov   2012")
    print(" ")


def analyze_geometry_block(f):
    pline = ""
    fido_per_cm = []
    for subline in f:
        if "T" not in subline:
            if "/" not in subline:
                if "=" not in subline:
                    pline += " " + subline.replace("\n", "").strip()
                else:
                    fine_check(pline)
                    mesh_boundary_check(pline)
                    # the pulling of FM per CM is actually done here
                    fido_per_cm = fido_grab(pline, fido_per_cm)
                    pline = ""
                    subline = subline.replace("\n", "").strip()
                    pline = "\n" + subline
        else:
            break
    return fido_per_cm


def print_fido_string(fido_per_cm):
    print("Printing FIDO strings for each coarse mesh"
          " (output may be truncated):")
    print("-"*80)
    string = ""
    for index, item in enumerate(fido_per_cm):
        if len(item) < 70:
            string += str(index + 1) + "-->" + item+"   "
        else:
            if len(string) > 0:
                print(textwrap.fill(string, width=80))
            string = str(index + 1) + "-->" + item[:70]+" ..."
            print(string)
            string = ""


def write_cell_card():
    print("Printing cell card")

    print("c MESH2MCNP MCNP Input Generator, Date generated: ",
          datetime.date.today(), file=output_fh)
    print("c input deck generated from PENTRAN input:",
          input_file, file=output_fh)
    ### cell cards
    print("c\nc cell cards -----\n" +
          "c description: each material\n" +
          "c   is assigned to a matching cell number\n" +
          "c   and matching universe (same as matid)\n" +
          "c   inside of the global rpp", file=output_fh)
    matlCell = produce_matl_cell()
    for cell in matlCell:
            print(cell, file=output_fh)

    print("c\nc cells 10001 to", str(10000+len(fido_per_cm)),
          " describe a cartesian f.m. lattice in each coarse mesh\nc",
          file=output_fh)


def write_cell_lattice(maxmat, fido_per_cm):
    ### This section pulls out the full material fine mesh per coarse mesh
    ### in fm_matl_per_cm
    fm_matl_per_cm = []
    for index, fido_string in enumerate(fido_per_cm):

        cellString = str(index + 1 + 10000).ljust(5) + "  0      -" + \
            str(index+1)+"\n"

        cellString += " "*9 + "lat=1" + " "*24 + "u=" + \
            str(index + 1 + 10000).ljust(5) + " imp:n=1 $ lattice\n"

        cellString += "         fill=0:"

        print(cellString +
              str(fm_attributes['xFm'][index]-1).ljust(4)+"0:" +
              str(fm_attributes['yFm'][index]-1).ljust(4)+"0:" +
              str(fm_attributes['zFm'][index]-1), file=output_fh)
        mcnp_string, materialString = process_pfm_string(fido_string)
        fm_matl_per_cm.append(materialString)

        print(mcnp_string, file=output_fh)

    print("c\nc cells 50001 to", str(50000+len(fido_per_cm)),
          " place the coarse meshes in universe 0, the real universe\nc",
          file=output_fh)

    for index in range(0, len(fido_per_cm)):  # could have used maxmat
        cellString = str(index + 1 + 50000).ljust(5) + "  0       -"\
            + str(index + 1 + 10000).ljust(5)\
            + " fill=" + str(index + 1 + 10000).ljust(5) + " imp:n=1"
        print(cellString, file=output_fh)

    print("c\nc a global rpp is defined +0.00000001",
          "from the coarse mesh boundaries\nc",
          " and also needs a non-zero importance assigned",
          file=output_fh)
    string = "99998  0  99990 99991 -99992 99993 -99994 99995 -99996 imp:n=1"

    print("c 99999 cell is the void region outside the problem boundaries",
          file=output_fh)
    print(string, file=output_fh)

    voidString = "99999  0 -99991:99992:-99993:99994:-99995:99996 imp:n=0"
    print(voidString, file=output_fh)

    return fm_matl_per_cm


def write_reference_file(fm_matl_per_cm, cm_attributes):

# note that coarse mesh indices are zero-based (c/python-style)
# dependencies of the loop below:
#   fm_matl_per_cm, x/y/zMinCm, ix/jy/kzfL

    reference_fh = open('prb.ref', 'w')
    xMinCm = cm_attributes['xMinCm']
    yMinCm = cm_attributes['yMinCm']
    zMinCm = cm_attributes['zMinCm']
    xFmDimL = cm_attributes['xFmDimL']
    yFmDimL = cm_attributes['xFmDimL']
    zFmDimL = cm_attributes['xFmDimL']
    kzfL = cm_attributes['kzfL']
    jyfL = cm_attributes['jyfL']
    ixfL = cm_attributes['ixfL']

    print('  Group       Coarse      Material   ',
          'x-cm        y-cm        z-cm        Phi0',
          '(needs to come from MCNP!)', file=reference_fh)

    for cm_idx in range(len(cm_attributes['fm_per_cm'])):

        matlPerCm = fm_matl_per_cm[cm_idx].split()

        beginTriplet = \
            [Decimal(xMinCm[cm_idx])+xFmDimL[cm_idx]/Decimal(2.),
             Decimal(yMinCm[cm_idx])+yFmDimL[cm_idx]/Decimal(2.),
             Decimal(zMinCm[cm_idx])+zFmDimL[cm_idx]/Decimal(2.)]

        count = 0
        for zdx in range(int(kzfL[cm_idx])):
            for ydx in range(int(jyfL[cm_idx])):
                for xdx in range(int(ixfL[cm_idx])):

                    masterTriplet = \
                        [beginTriplet[0]+xdx*xFmDimL[cm_idx],
                         beginTriplet[1]+ydx*yFmDimL[cm_idx],
                         beginTriplet[2]+zdx*zFmDimL[cm_idx]]

                    string = '  ' + '{:.4E}  ' * 6

                    print(string.format(1, cm_idx+1, float(matlPerCm[count]),
                          float(masterTriplet[0]),
                          float(masterTriplet[1]),
                          float(masterTriplet[2])), file=reference_fh)

                    count += 1

    reference_fh.close()


def write_surface_and_data(group_upper_boundaries):

    groups = instances["maxgrp"].value
    maxgcm = instances["maxgcm"].value
    maxmat = instances["maxmat"].value

    print("Printing surface card")
    ### surface card
    print("", file=output_fh)
    print("c surface cards -----", file=output_fh)
    rpp_list = produce_surface_rpp(cm_attributes)
    for rpp in rpp_list:
            print(rpp, file=output_fh)

    print("Printing data card")
    ### data card
    print("", file=output_fh)
    print("c data card -----", file=output_fh)
    print("mode n", file=output_fh)

    # material cards
    print("c material cards", file=output_fh)
    for i in range(1, maxmat+1):
        print("m" + str(i) + "  " + str(i) + "000.22m  1.0", file=output_fh)

    # multigroup option
    print("mgopt f " + str(groups), file=output_fh)

    # fmesh tally cards : option added nov.14.2012
    if group_upper_boundaries is not None:  # check if list is non-empty
        group_upper_boundaries.reverse()

    for i in range(maxgcm):
        xmin_by_cm = globals()['cm_attributes']['xMinCm']
        ymin_by_cm = globals()['cm_attributes']['yMinCm']
        zmin_by_cm = globals()['cm_attributes']['zMinCm']
        xmax_by_cm = globals()['cm_attributes']['xMaxCm']
        ymax_by_cm = globals()['cm_attributes']['yMaxCm']
        zmax_by_cm = globals()['cm_attributes']['zMaxCm']
        ixfl_by_cm = globals()['cm_attributes']['ixfL']
        jyfl_by_cm = globals()['cm_attributes']['jyfL']
        kzfl_by_cm = globals()['cm_attributes']['kzfL']
        # the actual printing for each fmesh
        print("fmesh" + str(i + 1)+"4:n  origin=", end="", file=output_fh)
        print(xmin_by_cm[i], ymin_by_cm[i], zmin_by_cm[i], file=output_fh)
        print('     ', 'imesh=', xmax_by_cm[i],
              'iints=', int(ixfl_by_cm[i]), file=output_fh)
        print('     ', 'jmesh=', ymax_by_cm[i],
              'jints=', int(jyfl_by_cm[i]), file=output_fh)
        print('     ', 'kmesh=', zmax_by_cm[i],
              'kints=', int(kzfl_by_cm[i]), file=output_fh)

        # the emesh card
        if not group_upper_boundaries:
            print("      emesh=", end="", file=output_fh)
            for i in range(groups):
                print("? ", end="", file=output_fh)
            print(" $ supply the",
                  input_file.split('.')[0] +
                  ".grp file and the '?' will go away!",
                  end="", file=output_fh)
            print(file=output_fh)
        else:
            group_string = "      emesh="
            for i in range(groups):
                group_string += str(group_upper_boundaries[i]) + " "

            print(textwrap.fill(group_string, width=80,
                  subsequent_indent='            '), file=output_fh)

        eint_string = "      eints="
        for i in range(groups):
            eint_string += "1 "

        print(textwrap.fill(eint_string, width=80,
              subsequent_indent='            '), file=output_fh)

    # f4 tally cards
    for i in range(1, maxgcm+1):
        index = i - 1
        print("c f"+str(i)+"4:n  ", end="", file=output_fh)
        tally_string = "( ("
        for m in range(1, maxmat+1):
            tally_string += " "+str(m)
        tally_string += " ) < "
        tally_string += str(10000 + i) + "[0:" + \
            str(fm_attributes['xFm'][index] - 1).ljust(4) + \
            "0:" + str(fm_attributes['yFm'][index] - 1).ljust(4) + \
            "0:" + str(fm_attributes['zFm'][index]-1) + "] < " + \
            str(50000 + i) + " )"
        print(tally_string, file=output_fh)


def insert_kcode():
    print("c if kcode or kcode.txt file exists,"
          " then its contents will be inserted here", file=output_fh)

    try:
        with open('kcode') as kfile:
            print("Appending contents of kcode...")
            for line in kfile:
                print(line, end="", file=output_fh)
    except IOError as e:  # noqa
        print('no \'kcode\' file detected, trying \'kcode.txt\'')
        try:
            with open('kcode.txt') as kfile:
                print("Appending contents of kcode.txt...")
                for line in kfile:
                    print(line, end="", file=output_fh)
        except IOError as e2:  # noqa
            print('')


def insert_group_boundaries(file):
    """ if the prb.grp file exists, read in the group information """

    group_upper_boundaries = []
    groups = instances["maxgrp"].value
    file_str = file.split('.')
    basename = file_str[0]
    try:
        with open(basename+'.grp') as gfile:
            print("Analyzing Group File for", groups, "groups...")
            next(gfile)
            next(gfile)
            next(gfile)
            next(gfile)
            gcnt = 0
            for line in gfile:
                #DBG print(line,end="")
                num, upperMeV_val = line.split()
                group_upper_boundaries.append(upperMeV_val)
                gcnt += 1
                if gcnt == groups:
                    break

        line_count = len(open(basename+'.grp', 'r').readlines())
        print("Group(.grp) file has", line_count, "lines")
    except IOError as e:  # noqa
        print("no grp file detected (optional)")

    return group_upper_boundaries

# main script begin

#lenArgv=len(sys.argv)
#if lenArgv == 1:
#        input_name=input("Enter PENTRAN file name: ")
#else:

input_file = sys.argv[1]
print(sys.argv[1])

#input_file=input("Enter PENTRAN file name: ")
#DBG input_file="uo2asm.pen"
output_file = input_file.split('.')[0]+'.mc'
print("Assigning name of output as:", output_file)

header()


input_fh = open(input_file, 'r')
#reference_file=open('prb.ref','w')

#GATHER INFORMATION

# initial search on fm, cm, geometry information in .pen file
for line in input_fh:

    # identify the maximum number of CM
    cm_check(input_fh, line)

    # identify the maximum number of grps, materials
    grp_matl_check(input_fh, line)

    # examine BLOCK II section
    if "BLOCK II(geometry)" in line:
        print("Analyzing Geometry Block...")
        fido_per_cm = analyze_geometry_block(input_fh)

    # identify boundary conditions
    bdy_check(input_fh, line)

print("\nPENTRAN Instances Summary")
print("-"*80, end="")
pretty_print_instances(instances)
print_fido_string(fido_per_cm)

fm_attributes = calc_fm_per_cm()
cm_attributes = \
    calc_boundary_per_cm(fm_attributes['xFm'],
                         fm_attributes['yFm'],
                         fm_attributes['zFm'])

### printing to output file
output_fh = open(output_file, 'w')

write_cell_card()
fm_matl_per_cm = write_cell_lattice(instances["maxmat"].value, fido_per_cm)

group_upper_boundaries = insert_group_boundaries(input_file)
write_surface_and_data(group_upper_boundaries)
insert_kcode()

print("\nPhysical Dimension Summary")
print("-"*80)
table_print('cm_attributes', 'xMinCm xMaxCm yMinCm yMaxCm zMinCm zMaxCm')
table_print('cm_attributes', 'xFmDimL yFmDimL zFmDimL ixfL jyfL kzfL')

write_reference_file(fm_matl_per_cm, cm_attributes)

output_fh.close()
