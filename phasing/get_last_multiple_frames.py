#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2018-Jun-12, latest modified on 2018-Jun-28

import argparse


def get_frame_num(inputfile): # get how many frame in inputfile (run.pdb), it doesnt used in this file
    framecount = 0
    with open (inputfile , 'r') as fr:
        for lines in fr.readlines():
            strs = lines.split()
            if strs[0] == 'END':
                framecount += 1
        fr.close()
    return framecount

def get_eachframe_lines(inputfile): # get how many lines for one frame
    frame_lines = 0
    with open (inputfile , 'r') as fr:
        for lines in fr.readlines():
            strs = lines.split()
            frame_lines += 1
            if strs[0] == 'END':
                break
        fr.close
    return frame_lines
        
def get_tail_frames(file, frame_lines, frame_number=20, avg_line_length=None):
    with open(file, 'r') as fr:
        if not avg_line_length: # Find average line length for next step
            fr.seek(0, 2)
            fr.seek(fr.tell() - 3000)
            avg_line_length = int(3000 / len(fr.readlines())) + 10
        
        fr.seek(0, 2)     # this block is for calculating offset value for pointer
        end_pointer = fr.tell()
        tail_lines = frame_number * frame_lines
        offset = tail_lines * avg_line_length
        if offset > end_pointer:
            fr.seek(0, 0)
            
        location = fr.tell() # get location now
        while True: # this part is the core part for moving pointer, use location to record position
            location = location - offset
            if location < 0:
                location = 0
                fr.seek(0, 0)
                break
            fr.seek(location)
            current_lines = len(fr.readlines())
            if current_lines >= tail_lines:
                break

        '''while len(fr.readlines()) < tail_lines: # this doesn't work if file length fewer than offset
            location = fr.tell() - offset
            fr.seek(location)
            if fr.tell() - offset < 0:
                fr.seek(0, 0)
                break'''

        fr.seek(location)
        checkpoint = fr.tell()
        
        lines = fr.readlines()
        index = 1
        
        fw = open('save' + str(index) + '.pdb', 'w')
        while index <= frame_number:
            for i in range(-tail_lines,0):
                if lines[i].split()[0] != 'END':
                    fw.write(lines[i])
                else:
                    fw.write(lines[i])
                    fw.close()
                    index += 1
                    if index > frame_number: # kill save21.pdb
                        break
                    fw = open('save' + str(index) + '.pdb', 'w')
        fw.close()
    fr.close()
    
    
def main():
    parser = argparse.ArgumentParser(
        description="This script extracts last several frames from a set of pdb file")
    parser.add_argument("input", help="inputfile")
    parser.add_argument("frames", help="number of frames you want", type=int)
    args = parser.parse_args()
    inputfile = args.input
    frames = args.frames
    
    #inputfile = 'run.pdb'
    #frames = 20
    
    frame_lines = get_eachframe_lines(inputfile)
    get_tail_frames(inputfile, frame_lines, frames)

if __name__== '__main__':
    main()
