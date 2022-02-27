# trussSTLGenerator

This converts a truss design (in terms of a list of coordinates and connections that include the cross sectional area) into STL files (either ascii or binary) with a predefined polygon count on the cross section.

'''
    input:
    fPath: path of the input file
    fCoord: file name of the coordinate file
    fCon: file name of the connectivity file
    fOutL: file name of the STL
    iDiv: number of polygon to approximate the circular cross section, 4 would mean a square
    boolAscii: True - ascii, false - binary (smaller file size but not human readable)
'''