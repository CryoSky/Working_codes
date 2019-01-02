#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2018-May-24, modified from Carlos's code, latest modified on 2018-Jun-13

import prody
import sys


folder = sys.argv[1]

best = prody.parsePDB('%s/1.pdb' % folder)
runs = []
best.delCoordset(0)

for i in range(1,21):
    try:
        temp = prody.parsePDB('%s/%i.pdb' % (folder,i))
        best.addCoordset(temp.getCoords())
        runs += [i]
    except ValueError as e:
        print (e,i)
        pass
    except OSError as e:
        print (e,i)
        pass

selection = best.select('name CA')

prody.alignCoordsets(best)
s = prody.parsePDB('%s/19.pdb' % folder)
s.setBetas(prody.calcRMSF(best))
prody.writePDB('%s/19_RMSF.pdb' % folder, s)

