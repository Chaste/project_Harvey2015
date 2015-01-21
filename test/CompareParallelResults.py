"""Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

#!/usr/bin/python

from itertools import izip

import math
import numpy
import sys
import matplotlib.pyplot as pp

# Check command line arguments
if (len(sys.argv) != 4):
  err = "Expected 3 command line arguments of SPACE_DIM and 2 .viznode files to compare, got " + str(len(sys.argv)-1);
  print err;
  sys.exit();

file_1_path = sys.argv[2];
file_2_path = sys.argv[3];

SPACE_DIM = int(sys.argv[1]);

DBL_EPS = numpy.finfo(numpy.float64).eps

file_1 = open(file_1_path, 'r');
file_2 = open(file_2_path, 'r');

num_lines_1 = sum(1 for line in file_1);
num_lines_2 = sum(1 for line in file_2);

if (num_lines_1 != num_lines_2):
  err = "Files differ in the number of lines, File 1 has " + str(num_lines_1) + " File 2 has " + str(num_lines_2);
  print err;
  sys.exit();

# Open files simultaneously
diff = 0.0;
test_points = 0;
time_difs = []

with open(file_1_path) as file_1:
  with open(file_2_path) as file_2: 
    for line_1, line_2 in izip(file_1, file_2):
      locations_1 = numpy.array(line_1.split());
      # Delete the time-stamp
      time = locations_1[0]
      locations_1 = locations_1[1:];
      # Reshape into list of node locations.
      locations_1 = numpy.reshape(locations_1, (-1, SPACE_DIM));
      locations_1 = numpy.sort(locations_1, 0);

      locations_2 = numpy.array(line_2.split());
      # Delete the time-stamp
      locations_2 = locations_2[1:];
      # Reshape into list of node locations.
      locations_2 = numpy.reshape(locations_2, (-1, SPACE_DIM));
      locations_2 = numpy.sort(locations_2, 0);

      local_diff = 0.0
      for x,y in izip(locations_1, locations_2):
        for i, j in izip(x,y):
          local_diff = local_diff + math.fabs(float(i) - float(j))
      sigma = local_diff / len(locations_1)
      #time_difs.append(sigma)
      time_difs.append(sigma /DBL_EPS)

output = ""
for idx, diff in enumerate(time_difs):
  output += "({a}, {b}) ".format(a=idx, b=diff)
print output 

pp.xlabel('t (hours)')
pp.ylabel('sigma_n / DBL_EPSILON')
pp.plot(time_difs)
pp.show()
