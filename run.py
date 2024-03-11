#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from mod import *

input_file = 'INPUT.json'

e, v0, atol, rtol, method, z1, z2, step, n_rot, save_dir = ReadParams(input_file)

CreateDir(save_dir)

zs = np.arange(z1, z2, step)
print('start solving phase map')
zv = []
for i in tqdm(zs):
    zv.append(PhaseMap(n_rot, i, v0, e, atol, rtol, method))

files = SavePhaseMap(e, z1, z2, step, zv, save_dir)

print('start plotting phase maps')
Plot(e, z2, save_dir, files)
print('done')