#!/usr/bin/python3

import os
import re

# run command
def run(command):
    print(command)
    return os.popen('/bin/bash -c \'{0}\''.format(command)).readlines()

# fl = run('find . -name 0806.log')

# s = ''
# for fn in fl:
#     if fn is None:
#         continue
#     fn = fn[:-1]
#     with open(fn, 'r') as f:
#         try:
#             while 1:
#                 line = next(f)
#                 if 'total score' in line:
#                     case = fn.split('/')[2]
#                     runtime = float(re.findall('\[(.+?)\]', line)[0].strip())
#                     score = float(re.findall('= (.+?)\n', line)[0])
#                     s += '{}:\t\t{:4.4f}    {:.1f}\n'.format(case, runtime, score)
#         except StopIteration:
#             pass

# print(s)
t1 = 0
t2 = 0
while 1:
    a = float(input())
    if a == 0:
        break
    b = int(input())
    print(a/b)
    t1 += a
    t2 += b
print(t1/t2)